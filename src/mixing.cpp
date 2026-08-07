#include "mixing.h"
#include "EOSlist.h"

//#include <gsl/gsl_math.h>
//#include <limits>

namespace Mixing {
  static const double EPS = 1e-12;

  // Clamp tiny negatives to zero, but fail for real negatives
  inline bool enforce_non_negative(double &v)
  {
    if (v < -EPS) return false;
    if (v < 0.0) v = 0.0;
    return true;
  }

  static bool g_upper_mantle_has_free_sio2 = false;

  bool upper_mantle_has_free_sio2()
  {
    return g_upper_mantle_has_free_sio2;
  }

  // Simple normalization helper
  //static void normalize(vector<double> &x)
  //{
  //  double sum = 0.0;
  //  for (size_t i = 0; i < x.size(); ++i)
  //    sum += x[i];

  //  if (sum <= 0.0) return;

  //  for (size_t i = 0; i < x.size(); ++i)
  //    x[i] /= sum;
  //}

  // -------------------- General ideal mixture density --------------------
  //
  // Ideal volume mixing:
  //   v_i  = M_i / rho_i     (cm^3 / mol)
  //   v_mix = Σ x_i v_i
  //   M_mix = Σ x_i M_i
  //   rho_mix = M_mix / v_mix
  //
  double density_ideal_mixture(double P_cgs, double T,
                               const vector<EOS*> &components,
                               const vector<double> &x_in, double rho_guess)
  {
    const size_t n = components.size();
    if (n == 0 || x_in.size() != n)
      return numeric_limits<double>::quiet_NaN();

    vector<double> x = x_in;
    //normalize(x);

    // sanitize guess (avoid absurd values propagating into component solvers)
    double guess = rho_guess;
    
    double V_mix = 0.0;   // cm^3/mol

    for (size_t i = 0; i < n; ++i)
    {
      if (x[i] <= 1e-10) continue;
      EOS *phase = components[i];
      const double Mi = phase->getmmol();        // g/mol
      if (!gsl_finite(Mi) || Mi <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      const double rho_i = phase->density(P_cgs, T, guess); // g/cm^3
      
      if (!gsl_finite(rho_i) || rho_i <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      V_mix += x[i] / rho_i;
    }

    if (V_mix <= 0.0)
      return numeric_limits<double>::quiet_NaN();

    return 1 / V_mix; // g/cm^3
  }

    // Overload: keep per-component rho guesses (same ordering as `components`)
  double density_ideal_mixture(double P_cgs, double T,
                               const vector<EOS*> &components,
                               const vector<double> &x_in,
                               vector<double> &rho_guess_phase,
                               double rho_guess_init)
  {
    const size_t n = components.size();
    if (n == 0 || x_in.size() != n)
      return numeric_limits<double>::quiet_NaN();

    vector<double> x = x_in;
    //normalize(x);

    if (!gsl_finite(rho_guess_init) || rho_guess_init <= 0.0)
      rho_guess_init = 1.0;

    // Initialize / sanitize guess array
    if (rho_guess_phase.size() != n)
      rho_guess_phase.assign(n, rho_guess_init);
    else
      for (size_t i = 0; i < n; ++i)
        if (!gsl_finite(rho_guess_phase[i]) || rho_guess_phase[i] <= 0.0)
          rho_guess_phase[i] = rho_guess_init;

    double V_mix = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
      if (x[i] <= 1e-10) continue;
      EOS *phase = components[i];

      const double rho_i = phase->density(P_cgs, T, rho_guess_phase[i]);
      if (!gsl_finite(rho_i) || rho_i <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      rho_guess_phase[i] = rho_i;   // <-- save for next call
      V_mix += x[i] / rho_i;
    }

    if (V_mix <= 0.0)
      return numeric_limits<double>::quiet_NaN();

    return 1.0 / V_mix;
  }

  // -------------------- General ideal mixture dT/dP|S --------------------
  //
  // Parallel addition of component gradients:
  //
  //   (dT/dP)_mix^{-1} = Σ x'_i / (dT/dP)_i
  //
  // where x'_i are the molar fractions renormalized over the
  // subset of components that actually have a thermal model
  // (thermal_type != 0).
  //
  double dTdP_S_ideal_mixture(double P_cgs, double T,
                              const vector<EOS*> &components,
                              const vector<double> &x_in, double rho_guess)
  {
    const size_t n = components.size();
    if (n == 0 || x_in.size() != n)
      return numeric_limits<double>::quiet_NaN();

    vector<double> x = x_in;
    //normalize(x);

    double inv_grad_sum = 0.0;  // Σ x_i / grad_i using original x_i
    double x_used_sum   = 0.0;  // Σ x_i over phases with a usable gradient
    int    n_thermal    = 0;    // number of phases contributing to gradient

    double guess = rho_guess;
    for (size_t i = 0; i < n; ++i)
    {
      if (x[i] <= 1e-10) continue;
      EOS *phase = components[i];

      // Skip EOS that have no thermal parameters: thermal_type == 0
      const int thermal_type = phase->getthermal();
      if (thermal_type == 0)
        continue;

      // We need a reasonable rho_guess for dTdP_S. Use pure-phase density.
      double rho_i = phase->density(P_cgs, T, guess);
      if (!gsl_finite(rho_i) || rho_i <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      const double grad_i = phase->dTdP_S(P_cgs, T, rho_i); // K/GPa

      // For EOS that *should* have a gradient, treat non-finite/zero as error.
      if (!gsl_finite(grad_i) || grad_i == 0.0)
        return numeric_limits<double>::quiet_NaN();

      inv_grad_sum += x[i] / grad_i;  // units: GPa/K
      x_used_sum   += x[i];
      ++n_thermal;
    }

    // If every component was thermally inert (thermal_type == 0),
    // keep the usual convention: isothermal (dT/dP = 0).
    if (n_thermal == 0)
      return 0.0;

    if (inv_grad_sum <= 0.0)
      return numeric_limits<double>::quiet_NaN();

    // Renormalization over thermal components:
    //
    //   x'_i = x_i / x_used_sum
    //   Σ x'_i / grad_i = (1 / x_used_sum) Σ x_i / grad_i
    //   => (dT/dP)_mix = x_used_sum / Σ (x_i / grad_i)
    //
    const double grad_mix = x_used_sum / inv_grad_sum;  // K/GPa
    return grad_mix;
  }

  // -------------------- General ideal mixture dT/dP|S (full rigorous version) --------------------
  //
  // More rigorous formulation accounting for volume expansion:
  //
  //   (dT/dP)_mix = (Σ x_i/m_i (∂V/∂T)_P,i) / (Σ x_i/m_i (∂V/∂T)_P,i / (dT/dP)_i)
  //
  // where (∂V/∂T)_P,i = (m_i/rho_i^2) * (pPpT_rho,i / pPprho_T,i)
  //
  // This accounts for the thermal expansion of each component.
  //
  double dTdP_S_ideal_mixture_full(double P_cgs, double T,
                                    const vector<EOS*> &components,
                                    const vector<double> &x_in, double rho_guess)
  {
    const size_t n = components.size();
    if (n == 0 || x_in.size() != n)
      return numeric_limits<double>::quiet_NaN();

    vector<double> x = x_in;
    //normalize(x);

    double P_GPa = P_cgs / 1E10;
    double numerator = 0.0;    // Σ x_i (∂V/∂T)_P,i
    double denominator = 0.0;  // Σ x_i (∂V/∂T)_P,i / (dT/dP)_i
    int    n_thermal = 0;      // number of phases contributing to gradient

    double guess = rho_guess;

    for (size_t i = 0; i < n; ++i)
    {
      if (x[i] <= 1e-10) continue;
      EOS *phase = components[i];

      // Skip EOS that have no thermal parameters: thermal_type == 0
      const int thermal_type = phase->getthermal();
      if (thermal_type == 0)
        continue;

      // Get component density
      double rho_i = phase->density(P_cgs, T, guess);
      if (!gsl_finite(rho_i) || rho_i <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      // Get molar mass
      const double m_i = phase->getmmol(); // g/mol
      if (!gsl_finite(m_i) || m_i <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      // Calculate (∂V/∂T)_P,i using the dVdT_P function
      // For 4-column tables (eqntype==7 && thermal_type==10), this directly uses table data
      // For other types, it uses: (m_i/rho_i^2) * (pPpT_rho,i / pPprho_T,i)
      const double dVdT_P_i = phase->dVdT_P(P_GPa, T, rho_i);

      if (!gsl_finite(dVdT_P_i))
        return numeric_limits<double>::quiet_NaN();

      // Get component adiabatic gradient
      const double grad_i = phase->dTdP_S(P_cgs, T, rho_i); // K/GPa

      // For EOS that *should* have a gradient, treat non-finite/zero/negative as error.
      // Negative grad_i would cause denominator to become negative, leading to large negative dTdP
      if (!gsl_finite(grad_i) || grad_i <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      // Accumulate terms
      numerator += x[i]/m_i * dVdT_P_i;
      denominator += x[i]/m_i * dVdT_P_i / grad_i;
      ++n_thermal;
    }

    // If every component was thermally inert (thermal_type == 0),
    // keep the usual convention: isothermal (dT/dP = 0).
    if (n_thermal == 0)
      return 0.0;

    if (denominator <= 0.0 || !gsl_finite(numerator) || !gsl_finite(denominator))
      return numeric_limits<double>::quiet_NaN();

    // Calculate final gradient: (dT/dP)_mix = numerator / denominator
    const double grad_mix = numerator / denominator;  // K/GPa
    return grad_mix;
  }

    // Overload: keep per-component rho guesses (same ordering as `components`)
  double dTdP_S_ideal_mixture_full(double P_cgs, double T,
                                   const vector<EOS*> &components,
                                   const vector<double> &x_in,
                                   vector<double> &rho_guess_phase,
                                   double rho_guess_init,
                                   bool use_cached_rho = false)
  {
    const size_t n = components.size();
    if (n == 0 || x_in.size() != n)
      return numeric_limits<double>::quiet_NaN();

    vector<double> x = x_in;
    //normalize(x);

    if (!gsl_finite(rho_guess_init) || rho_guess_init <= 0.0)
      rho_guess_init = 1.0;

    if (rho_guess_phase.size() != n)
      rho_guess_phase.assign(n, rho_guess_init);
    else
      for (size_t i = 0; i < n; ++i)
        if (!gsl_finite(rho_guess_phase[i]) || rho_guess_phase[i] <= 0.0)
          rho_guess_phase[i] = rho_guess_init;

    const double P_GPa = P_cgs / 1E10;

    double numerator   = 0.0;
    double denominator = 0.0;
    int n_thermal = 0;

    for (size_t i = 0; i < n; ++i)
    {
      if (x[i] <= 1e-10) continue;
      EOS *phase = components[i];

      const int thermal_type = phase->getthermal();
      if (thermal_type == 0)
        continue;

      double rho_i;
      if (use_cached_rho && gsl_finite(rho_guess_phase[i]) && rho_guess_phase[i] > 0.0)
      {
        rho_i = rho_guess_phase[i];   // reuse density from preceding density_NAME call
      }
      else
      {
        rho_i = phase->density(P_cgs, T, rho_guess_phase[i]);
        if (!gsl_finite(rho_i) || rho_i <= 0.0)
          return numeric_limits<double>::quiet_NaN();
        rho_guess_phase[i] = rho_i;
      }

      const double m_i = phase->getmmol();
      if (!gsl_finite(m_i) || m_i <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      // IMPORTANT: use the rho-taking overload so dVdT_P doesn't call density() again
      const double dVdT_P_i = phase->dVdT_P(P_GPa, T, rho_i);
      if (!gsl_finite(dVdT_P_i))
        return numeric_limits<double>::quiet_NaN();

      const double grad_i = phase->dTdP_S(P_cgs, T, rho_i); // K/GPa
      if (!gsl_finite(grad_i) || grad_i <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      numerator   += x[i] / m_i * dVdT_P_i;
      denominator += x[i] / m_i * dVdT_P_i / grad_i;
      ++n_thermal;
    }

    if (n_thermal == 0)
      return 0.0;

    if (denominator <= 0.0 || !gsl_finite(numerator) || !gsl_finite(denominator))
      return numeric_limits<double>::quiet_NaN();

    return numerator / denominator;
  }

  // --- Generic wrapper generator for a given mixture NAME ---
  //
  // Expects:
  //   vector<EOS*> comps_NAME;
  //   vector<double> x_NAME;
  //
  // and generates:
  //   density_NAME(P_cgs, T, rho_guess)
  //   dTdP_S_NAME(P_cgs, T, rho_guess)
  //   dTdP_NAME(P_cgs, T, rho_guess)
  //
  #define DEFINE_IDEAL_MIX_WRAPPERS(NAME)                                       \
  static vector<double> rho_guess_phase_##NAME;                                 \
  static double last_P_##NAME = -1;                                             \
  static double last_T_##NAME = -1;                                             \
                                                                               \
  double density_##NAME(double P_cgs, double T, double rho_guess)               \
  {                                                                             \
    double result = density_ideal_mixture(P_cgs, T, comps_##NAME, x_##NAME,     \
                                         rho_guess_phase_##NAME, rho_guess);   \
    if (gsl_finite(result)) {                                                   \
      last_P_##NAME = P_cgs;                                                    \
      last_T_##NAME = T;                                                        \
    }                                                                           \
    return result;                                                              \
  }                                                                             \
                                                                               \
  double dTdP_S_##NAME(double P_cgs, double T, double &rho_guess)               \
  {                                                                             \
    bool cached = (P_cgs == last_P_##NAME && T == last_T_##NAME);               \
    return dTdP_S_ideal_mixture_full(P_cgs, T, comps_##NAME, x_##NAME,          \
                                    rho_guess_phase_##NAME, rho_guess, cached);\
  }                                                                             \
                                                                               \
  double dTdP_##NAME(double P_cgs, double T, double &rho_guess)                 \
  {                                                                             \
    bool cached = (P_cgs == last_P_##NAME && T == last_T_##NAME);               \
    const double grad_GPa = dTdP_S_ideal_mixture_full(P_cgs, T,                 \
                                comps_##NAME, x_##NAME, rho_guess_phase_##NAME, \
                                rho_guess, cached);                            \
    if (!gsl_finite(grad_GPa))                                                  \
      return numeric_limits<double>::quiet_NaN();                               \
    return grad_GPa;                                                           \
  }

  bool mantle_wFeO_from_MgNumber(
  double CaMg, double SiMg, double AlMg,
  double Mg_number,
  double &mantle_wFeO_out,
  std::string &note_or_error)
  {
  note_or_error.clear();
  std::stringstream msg;

  if (!gsl_finite(Mg_number) || Mg_number <= 0.0 || Mg_number > 1.0) {
    note_or_error = "mantle_Mg_number must be in (0,1].";
    return false;
  }

  // Same molar-mass bookkeeping used by compute_mode9_core_mantle()
  const double mu_MgO   = mMg + mO;
  const double mu_SiO2  = mSi + 2.0*mO;
  const double mu_CaO   = mCa + mO;
  const double mu_AlO15 = mAl + 1.5*mO;
  const double mu_FeO   = mFe + mO;

  // Fe-free mantle mass per 1 mol Mg (same as mode9 core/mantle solver)
  const double mu_bar = mu_MgO + SiMg*mu_SiO2 + CaMg*mu_CaO + AlMg*mu_AlO15;
  if (!gsl_finite(mu_bar) || mu_bar <= 0.0) {
    msg << "Invalid mu_bar=" << mu_bar << " from ratios.";
    note_or_error = msg.str();
    return false;
  }

  const double k = mu_bar / mu_FeO;

  // Mg# = Mg/(Mg+Fe)  =>  Fe/Mg = (1-Mg#)/Mg#
  const double FeMg_mantle = (1.0 - Mg_number) / Mg_number;

  // wFeO = FeMg / (FeMg + k)  (algebraic inverse of compute_FeMg_mantle())
  mantle_wFeO_out = FeMg_mantle / (FeMg_mantle + k);

  if (!gsl_finite(mantle_wFeO_out) || mantle_wFeO_out < 0.0 || mantle_wFeO_out >= 1.0) {
    msg << "No physical mantle_wFeO from mantle Mg#=" << Mg_number
        << " (computed wFeO=" << mantle_wFeO_out << ").";
    note_or_error = msg.str();
    return false;
  }

  msg << "Converted mantle Mg#=" << Mg_number
      << " -> mantle_wFeO=" << mantle_wFeO_out
      << " (mantle Fe/Mg=" << FeMg_mantle << ").";
  note_or_error = msg.str();
  return true;
  }

  bool compute_mode9_core_mantle(
    double CaMg, double SiMg, double AlMg,
    double &FeMg_bulk_io,       // <0 => auto
    double &mantle_wFeO_io,     // <0 => auto
    double &RCMF_io,            // <0 => auto
    double wt_fract_S_core,
    double wt_fract_O_core,
    double wt_fract_H_core,
    double wt_fract_C_core,
    double wt_fract_Si_core,
    double wt_fract_Ni_core,
    double &FeMg_mantle_out,
    std::string &note_or_error)
  {
  note_or_error.clear();
  std::stringstream msg;

  auto is_set = [&](double v){ return v >= 0.0; };

  int nset = (int)is_set(FeMg_bulk_io) + (int)is_set(mantle_wFeO_io) + (int)is_set(RCMF_io);
  if (nset < 2) { note_or_error = "Need any 2 of {FeMg_bulk, mantle_wFeO, RCMF}. Set missing one to -1."; return false; }

  // ---- validate fixed core fractions ----
  if (wt_fract_S_core < 0.0)  wt_fract_S_core = 0.0;
  if (wt_fract_Ni_core < 0.0) wt_fract_Ni_core = 0.0;
  if (wt_fract_O_core < 0.0) wt_fract_O_core = 0.0;
  if (wt_fract_H_core < 0.0) wt_fract_H_core = 0.0;
  if (wt_fract_C_core < 0.0) wt_fract_C_core = 0.0;
  if (wt_fract_Si_core < 0.0) wt_fract_Si_core = 0.0;

  const double w_nonFe = wt_fract_S_core +wt_fract_O_core +wt_fract_H_core +wt_fract_C_core +wt_fract_Si_core +wt_fract_Ni_core;
  if (w_nonFe >= 1.0 - EPS) {
    msg << "Core non-Fe mass fraction =" << w_nonFe << " must be < 1.";
    note_or_error = msg.str();
    return false;
  }

  // ---- molar masses (g/mol) ----
  const double mu_MgO   = mMg + mO;
  const double mu_SiO2  = mSi + 2.0*mO;
  const double mu_CaO   = mCa + mO;
  const double mu_AlO15 = mAl + 1.5*mO;
  const double mu_FeO   = mFe + mO;
  const double mu_Fe    = mFe;

  // ---- mu_bar (Fe-free mantle mass per 1 mol Mg) ----
  const double mu_bar = mu_MgO + SiMg*mu_SiO2 + CaMg*mu_CaO + AlMg*mu_AlO15;
  if (!gsl_finite(mu_bar) || mu_bar <= 0.0) {
    msg << "Invalid mu_bar=" << mu_bar << " from ratios.";
    note_or_error = msg.str();
    return false;
  }

  const double A = mu_Fe / (1.0 - w_nonFe);   // core mass per mol Fe, accounting for non-Fe mass fraction
  const double k = mu_bar / mu_FeO;

  auto compute_FeMg_mantle = [&](double wFeO)->double{
    if (wFeO <= 0.0) return 0.0;
    return (wFeO / (1.0 - wFeO)) * k;
  };

  auto compute_R_from_F_and_w = [&](double FeMg_bulk, double wFeO, double &FeMg_mantle)->double{
    FeMg_mantle = compute_FeMg_mantle(wFeO);
    const double nFe_core = FeMg_bulk - FeMg_mantle;
    if (nFe_core <= 0.0) return 0.0;
    const double Mcore  = nFe_core * A;
    const double Mmant  = mu_bar / (1.0 - wFeO);
    return Mcore / (Mcore + Mmant);
  };

  // ---- solve missing parameter ----
  if (!is_set(FeMg_bulk_io)) {
    // given RCMF and wFeO -> solve FeMg_bulk
    const double wFeO = mantle_wFeO_io;
    const double R    = RCMF_io;

    if (wFeO < 0.0 || wFeO >= 1.0) { note_or_error="mantle_wFeO must be in [0,1)."; return false; }
    if (R <= 0.0 || R >= 1.0)      { note_or_error="RCMF must be in (0,1)."; return false; }

    const double FeMg_mantle = compute_FeMg_mantle(wFeO);
    const double Mmant = mu_bar / (1.0 - wFeO);
    const double nFe_core = (R * Mmant) / (A * (1.0 - R));

    FeMg_bulk_io = nFe_core + FeMg_mantle;

  } else if (!is_set(mantle_wFeO_io)) {
    // given RCMF and FeMg_bulk -> solve wFeO (analytic)
    const double FeMg_bulk = FeMg_bulk_io;
    const double R         = RCMF_io;

    if (FeMg_bulk <= 0.0)     { note_or_error="FeMg_bulk must be > 0."; return false; }
    if (R <= 0.0 || R >= 1.0) { note_or_error="RCMF must be in (0,1)."; return false; }

    const double C = (R * mu_bar) / (A * (1.0 - R));
    const double wFeO = (FeMg_bulk - C) / (FeMg_bulk + k);

    if (!gsl_finite(wFeO) || wFeO < 0.0 || wFeO >= 1.0) {
      msg << "No physical mantle_wFeO solution from FeMg_bulk=" << FeMg_bulk
          << " and RCMF=" << R << " (computed wFeO=" << wFeO << ").";
      note_or_error = msg.str();
      return false;
    }

    mantle_wFeO_io = wFeO;

  } else if (!is_set(RCMF_io)) {
    // given FeMg_bulk and wFeO -> solve RCMF
    const double FeMg_bulk = FeMg_bulk_io;
    const double wFeO      = mantle_wFeO_io;

    if (FeMg_bulk <= 0.0)        { note_or_error="FeMg_bulk must be > 0."; return false; }
    if (wFeO < 0.0 || wFeO >= 1.0) { note_or_error="mantle_wFeO must be in [0,1)."; return false; }

    double FeMg_mantle_tmp = 0.0;
    RCMF_io = compute_R_from_F_and_w(FeMg_bulk, wFeO, FeMg_mantle_tmp);

  } else {
    // all 3 provided: just consistency-check (optional)
    double FeMg_mantle_tmp = 0.0;
    const double Rcalc = compute_R_from_F_and_w(FeMg_bulk_io, mantle_wFeO_io, FeMg_mantle_tmp);
    if (fabs(Rcalc - RCMF_io) > 5e-3) {
      msg << "Warning: Provided (FeMg_bulk, mantle_wFeO, RCMF) are inconsistent. "
          << "Computed RCMF=" << Rcalc << " from FeMg_bulk and mantle_wFeO. ";
    }
  }

  // ---- final outputs ----
  FeMg_mantle_out = compute_FeMg_mantle(mantle_wFeO_io);

  const double nFe_core = FeMg_bulk_io - FeMg_mantle_out;
  if (nFe_core < -EPS) {
    msg << "Derived mantle Fe/Mg=" << FeMg_mantle_out
        << " exceeds bulk Fe/Mg=" << FeMg_bulk_io << ".";
    note_or_error = msg.str();
    return false;
  }

  msg << "Solved mode9 core/mantle: FeMg_bulk=" << FeMg_bulk_io
      << ", mantle_wFeO=" << mantle_wFeO_io
      << ", RCMF=" << RCMF_io
      << ", FeMg_mantle=" << FeMg_mantle_out
      << ", core non-Fe=" << w_nonFe << ".";
  note_or_error = msg.str();
  return true;
  }


  
/*
 * upper_out (size 13):
 *   0: Fo          (Mg2SiO4)
 *   1: Fa          (Fe2SiO4)
 *   2: En          (Mg2Si2O6)
 *   3: Fs          (Fe2Si2O6)
 *   4: St          (SiO2)
 *   5: MgO         (periclase)
 *   6: FeO         (wustite)
 *   7: Di          (CaMgSi2O6)
 *   8: Hd          (CaFeSi2O6)
 *   9: Py          (Mg3Al2Si3O12)
 *  10: Alm         (Fe3Al2Si3O12)
 *  11: Spinel       (MgAl2O4)
 *  12: Hercynite     (FeAl2O4)
 */
bool compute_upper_mantle_fractions(double CaMg,
                                    double SiMg,
                                    double AlMg,
                                    double FeMg,
                                    vector<double> &upper_out)
{
    upper_out.assign(13, 0.0);

    if (CaMg < 0.0 || AlMg < 0.0 || FeMg < 0.0 || SiMg <= 0.0)
        return false;

    const double R_Ca = CaMg;
    const double R_Si = SiMg;
    const double R_Al = AlMg;
    const double R_Fe = FeMg;

    const double X_Mg = 1.0 / (1.0 + R_Fe);
    const double X_Fe = 1.0 - X_Mg;

    // Fixed Ca/Al silicates (always formed)
    const double n_Cpx = R_Ca;        // Ca(Mg,Fe)Si2O6

    // Budgets remaining after Ca/Al silicates:
    //   Cpx uses: 1 (Mg+Fe) and 2 Si
    //   Decide if Al hosted in Garnets or Spinels
    //   Gt  uses: 3 (Mg+Fe) and 3 Si
    //   Sp uses: 1 Mg 
    const double n_Gt_if_garnet = 0.5 * R_Al;

    const double rem_Si_if_garnet =
        R_Si - (2.0*n_Cpx + 3.0*n_Gt_if_garnet);

    const double rem_cat_if_garnet =
        (1.0 + R_Fe) - (1.0*n_Cpx + 3.0*n_Gt_if_garnet);

    const bool use_spinel_for_al =
        (rem_Si_if_garnet  < 0.5*rem_cat_if_garnet - EPS) ||
        (rem_Si_if_garnet  < -EPS) ||
        (rem_cat_if_garnet < -EPS);

    // Al host:
    //   Garnet:  (Mg,Fe)3Al2Si3O12 -> 0.5 Al/Mg formula units
    //   Spinel:  (Mg,Fe)Al2O4      -> 0.5 Al/Mg formula units
    const double n_Gt = use_spinel_for_al ? 0.0 : 0.5 * R_Al;
    const double n_Sp = use_spinel_for_al ? 0.5 * R_Al : 0.0;

    // Fixed phases before Mg-Si balancing.
    const double Si_fixed  = 2.0*n_Cpx + 3.0*n_Gt;
    const double Cat_fixed = 1.0*n_Cpx + 3.0*n_Gt + 1.0*n_Sp;

    double rem_Si  = R_Si - Si_fixed;
    double rem_cat = (1.0 + R_Fe) - Cat_fixed;

    if (rem_Si  < -EPS) return false;
    if (rem_cat < -EPS) return false;
    if (rem_Si  < 0.0) rem_Si  = 0.0;
    if (rem_cat < 0.0) rem_cat = 0.0;

    // Phase pools, formula-unit moles.
    double n_Ol  = 0.0;  // (Mg,Fe)2SiO4
    double n_Opx = 0.0;  // (Mg,Fe)2Si2O6
    double n_St  = 0.0;  // SiO2
    double n_Ox  = 0.0;  // (Mg,Fe)O

    // Reduced-budget Mg-Si switch.
    if (rem_Si > rem_cat + EPS) {
        // High Si: no olivine. Use all cations in Opx; excess Si -> SiO2.
        n_Ol  = 0.0;
        n_Opx = 0.5 * rem_cat;
        n_St  = rem_Si - rem_cat;
        n_Ox  = 0.0;
    } else if (rem_Si < 0.5*rem_cat - EPS) {
        // Low Si / high Mg: no Opx. Use all Si in Ol; excess cations -> oxides.
        n_Opx = 0.0;
        n_Ol  = rem_Si;
        n_Ox  = rem_cat - 2.0*rem_Si;
        n_St  = 0.0;
    } else {
        // Nominal: Ol + Opx.
        n_Ol  = rem_cat - rem_Si;
        n_Opx = rem_Si - 0.5*rem_cat;
        n_St  = 0.0;
        n_Ox  = 0.0;
    }

    if (!enforce_non_negative(n_Ol))  return false;
    if (!enforce_non_negative(n_Opx)) return false;
    if (!enforce_non_negative(n_St))  return false;
    if (!enforce_non_negative(n_Ox))  return false;

    // Endmember moles.
    const double n_Fo  = X_Mg * n_Ol;
    const double n_Fa  = X_Fe * n_Ol;

    const double n_En  = X_Mg * n_Opx;
    const double n_Fs  = X_Fe * n_Opx;

    const double n_MgO = X_Mg * n_Ox;
    const double n_FeO = X_Fe * n_Ox;

    // Periclase/Wustite EOSs are stored as Mg4O4 / Fe4O4
    const double n_MgO_EOS = 0.25 * n_MgO;
    const double n_FeO_EOS = 0.25 * n_FeO;

    const double n_Di  = X_Mg * n_Cpx;
    const double n_Hd  = X_Fe * n_Cpx;

    const double n_Py  = X_Mg * n_Gt;
    const double n_Alm = X_Fe * n_Gt;

    const double n_Spinel    = X_Mg * n_Sp;
    const double n_Hercynite = X_Fe * n_Sp;

    const double n[13] = {
      n_Fo, n_Fa,
      n_En, n_Fs,
      n_St,
      n_MgO_EOS, n_FeO_EOS,
      n_Di, n_Hd,
      n_Py, n_Alm,
      n_Spinel, n_Hercynite
    };

    const double M[13] = {
      Forsterite->getmmol(),
      Fayalite->getmmol(),
      Enstatite->getmmol(),
      Ferrosilite->getmmol(),
      Stishovite->getmmol(), 
      Periclase->getmmol(),
      Wustite->getmmol(),
      Diopside->getmmol(),
      Hedenbergite->getmmol(),
      Pyrope->getmmol(),
      Almandine->getmmol(),
      Spinel->getmmol(),
      Hercynite->getmmol()
    };

    double mtot = 0.0;
    for (int i = 0; i < 13; ++i) mtot += n[i] * M[i];
    if (mtot <= 0.0) return false;

    for (int i = 0; i < 13; ++i) upper_out[i] = (n[i] * M[i]) / mtot;

    return true;
}

/*
 * middle_out (size 10):
 *   0: Wadsleyite (Mg)   (or Ringwoodite)
 *   1: Fe-Wadsleyite     (or Fe-Ringwoodite)
 *   2: Stishovite        (SiO2)
 *   3: Periclase         (MgO)
 *   4: Wustite           (FeO)
 *   5: Pyrope            (Mg3Al2Si3O12)
 *   6: Almandine         (Fe3Al2Si3O12)
 *   7: Grossular         (Ca3Al2Si3O12)
 *   8: Mg-Majorite       (Mg4Si4O12)
 *   9: Fe-Majorite*      (FeSiO3 placeholder via Fe_Akimotoite)
 *
 * The middle-mantle stoichiometry is solved in MgSiO3-equivalent units:
 *   majorite_eq = (Mg,Fe)SiO3
 *
 * but mapped to EOS formula units as:
 *   Mg-majorite EOS units  = Mg-majorite equivalent units / 4
 *   Fe-majorite EOS units  = Fe-majorite equivalent units
 *
 * Switch logic (after forming Ca/Al garnets):
 *   - Low Si:  Wds + (MgO/FeO), no majorite
 *   - Nominal: Wds + majorite_eq
 *   - High Si: majorite_eq + St, no Wds
 *
 * Assumption: as in the upper mantle, all Fe-Mg solid solutions in this region
 * share the bulk Mg-number X_Mg = 1/(1+Fe/Mg).
 */
bool compute_middle_mantle_fractions(double CaMg,
                                     double SiMg,
                                     double AlMg,
                                     double FeMg,
                                     vector<double> &middle_out)
{
    middle_out.assign(10, 0.0);

    if (CaMg < 0.0 || AlMg < 0.0 || FeMg < 0.0 || SiMg <= 0.0)
        return false;

    const double R_Ca = CaMg;
    const double R_Si = SiMg;
    const double R_Al = AlMg;
    const double R_Fe = FeMg;

    const double X_Mg = 1.0 / (1.0 + R_Fe);
    const double X_Fe = 1.0 - X_Mg;

    // Fixed Ca/Al garnets
    double n_Gross = R_Ca / 3.0;                // grossular
    double n_G     = 0.5 * R_Al - R_Ca / 3.0;   // (pyrope+almandine) pool

    if (!enforce_non_negative(n_Gross)) return false;
    if (!enforce_non_negative(n_G))     return false;

    // Remaining budgets after Ca/Al garnets:
    //   Si used = 3*(n_Gross + n_G) = 1.5*R_Al
    //   (Mg+Fe) cations used = 3*n_G
    double rem_Si  = R_Si - 3.0*(n_Gross + n_G);
    double rem_cat = (1.0 + R_Fe) - 3.0*n_G;

    if (rem_Si  < -EPS) return false;
    if (rem_cat < -EPS) return false;
    if (rem_Si  < 0.0) rem_Si  = 0.0;
    if (rem_cat < 0.0) rem_cat = 0.0;

    // Phase pools
    double n_W      = 0.0;  // total Wds/Rwd formula units: (Mg,Fe)2SiO4
    double n_Mj_eq  = 0.0;  // majorite reservoir in (Mg,Fe)SiO3-equivalent units
    double n_Ox     = 0.0;  // (Mg,Fe)O pool
    double n_St     = 0.0;  // SiO2

    // Reduced-budget nominal window:
    //   rem_cat/2 <= rem_Si <= rem_cat
    //
    // Stoichiometry in equivalent units:
    //   rem_Si  = n_W + n_Mj_eq + n_St
    //   rem_cat = 2*n_W + n_Mj_eq + n_Ox
    if (rem_Si < 0.5*rem_cat - EPS) {
        // Low Si: no majorite; Wds + oxides
        n_W     = rem_Si;
        n_Mj_eq = 0.0;
        n_Ox    = rem_cat - 2.0*n_W;
        n_St    = 0.0;
    } else if (rem_Si <= rem_cat + EPS) {
        // Nominal: Wds + majorite_eq
        n_W     = rem_cat - rem_Si;
        n_Mj_eq = 2.0*rem_Si - rem_cat;
        n_Ox    = 0.0;
        n_St    = 0.0;
    } else {
        // High Si: no Wds; majorite_eq + St
        n_W     = 0.0;
        n_Mj_eq = rem_cat;
        n_Ox    = 0.0;
        n_St    = rem_Si - rem_cat;
    }

    if (!enforce_non_negative(n_W))     return false;
    if (!enforce_non_negative(n_Mj_eq)) return false;
    if (!enforce_non_negative(n_Ox))    return false;
    if (!enforce_non_negative(n_St))    return false;

    // Split pools into endmembers using the bulk Mg#
    const double n_Wds      = X_Mg * n_W;
    const double n_FeWds    = X_Fe * n_W;
    const double n_MgO      = X_Mg * n_Ox;
    const double n_FeO      = X_Fe * n_Ox;
    const double n_MgO_EOS = 0.25 * n_MgO;
    const double n_FeO_EOS = 0.25 * n_FeO;  
    const double n_Py       = X_Mg * n_G;
    const double n_Alm      = X_Fe * n_G;

    // Majorite reservoir:
    //   Mg equivalent units -> Mg4Si4O12 EOS formula units by /4
    //   Fe equivalent units -> FeSiO3 EOS formula units directly
    const double n_MgMaj_eq   = X_Mg * n_Mj_eq;
    const double n_FeMaj_eq   = X_Fe * n_Mj_eq;
    const double n_MgMaj_EOS  = 0.25 * n_MgMaj_eq;
    const double n_FeMaj_EOS  = n_FeMaj_eq;

    // Convert formula-unit moles -> mass fractions
    const double n[10] = {
      n_Wds, n_FeWds, n_St, n_MgO_EOS, n_FeO_EOS,
      n_Py, n_Alm, n_Gross, n_MgMaj_EOS, n_FeMaj_EOS
    };

    const double M[10] = {
      Mg_Wadsleyite->getmmol(),
      Fe_Wadsleyite->getmmol(),
      Stishovite->getmmol(),
      Periclase->getmmol(),
      Wustite->getmmol(),
      Pyrope->getmmol(),
      Almandine->getmmol(),
      Grossular->getmmol(),
      Mg_Majorite->getmmol(),
      Fe_Akimotoite->getmmol()
    };

    double mtot = 0.0;
    for (int i = 0; i < 10; ++i) mtot += n[i] * M[i];
    if (mtot <= 0.0) return false;

    for (int i = 0; i < 10; ++i) middle_out[i] = (n[i] * M[i]) / mtot;

    return true;
}


/*
 * lower_out (size 7):
 *   0: Mg-Bridgmanite  (MgSiO3) (or Mg-PPv)
 *   1: Fe-Bridgmanite  (FeSiO3) (or Fe-PPv)
 *   2: Al-Bridgmanite  (AlAlO3)   [Al2O3 component]
 *   3: Stishovite     (SiO2)
 *   4: Periclase      (MgO)
 *   5: Wustite        (FeO)
 *   6: Ca-Perovskite  (CaSiO3)
 *
 * Switch logic:
 *   - Low Si is naturally handled by (MgO/FeO) pool.
 *   - High Si: cap Pv by available (Mg+Fe) cations; excess Si -> St.
 */
bool compute_lower_mantle_fractions(double CaMg,
                                    double SiMg,
                                    double AlMg,
                                    double FeMg,
                                    vector<double> &lower_out)
{
    lower_out.assign(7, 0.0);

    if (CaMg < 0.0 || AlMg < 0.0 || FeMg < 0.0 || SiMg <= 0.0)
        return false;

    const double R_Ca = CaMg;
    const double R_Si = SiMg;
    const double R_Al = AlMg;
    const double R_Fe = FeMg;

    const double X_Mg = 1.0 / (1.0 + R_Fe);
    const double X_Fe = 1.0 - X_Mg;

    double n_CaPv = R_Ca;
    double n_AlPv = 0.5 * R_Al;   // Al2O3 component in Pv/PPv (does not change Pv/Fp budgets here)

    // Requested Pv from Si (same as original scheme)
    double n_Pv_req = R_Si - R_Ca;
    if (!enforce_non_negative(n_Pv_req)) return false;

    // Total (Mg+Fe) cation capacity for Pv + oxides
    const double cap = 1.0 / X_Mg;  // = 1 + R_Fe

    double n_Pv = 0.0;
    double n_Fp = 0.0;  // (MgO+FeO) pool
    double n_St = 0.0;

    if (n_Pv_req > cap + EPS) {
        // High Si: cap Pv, dump excess Si into St
        n_Pv = cap;
        n_Fp = 0.0;
        n_St = n_Pv_req - cap;
    } else {
        n_Pv = n_Pv_req;
        n_Fp = cap - n_Pv;
        n_St = 0.0;
    }

    if (!enforce_non_negative(n_CaPv)) return false;
    if (!enforce_non_negative(n_AlPv)) return false;
    if (!enforce_non_negative(n_Pv))   return false;
    if (!enforce_non_negative(n_Fp))   return false;
    if (!enforce_non_negative(n_St))   return false;

    const double n_PvMg = X_Mg * n_Pv;
    const double n_PvFe = X_Fe * n_Pv;
    const double n_FpMg = X_Mg * n_Fp;
    const double n_FpFe = X_Fe * n_Fp;
    const double n_FpMg_EOS = 0.25 * n_FpMg;
    const double n_FpFe_EOS = 0.25 * n_FpFe;

    const double n[7] = { n_PvMg, n_PvFe, n_AlPv, n_St, n_FpMg_EOS, n_FpFe_EOS, n_CaPv };
    const double M[7] = {
      Mg_Bridgmanite->getmmol(),
      Fe_Bridgmanite->getmmol(),
      Al_Bridgmanite->getmmol(),
      Stishovite->getmmol(),
      Periclase->getmmol(),
      Wustite->getmmol(),
      Ca_Perovskite->getmmol()
    };

    double mtot = 0.0;
    for (int i = 0; i < 7; ++i) mtot += n[i] * M[i];
    if (mtot <= 0.0) return false;

    for (int i = 0; i < 7; ++i) lower_out[i] = (n[i] * M[i]) / mtot;

    return true;
}

/*
 * Mantle stoichiometry wrapper with minimal clamping + clear warnings:
 *  - clamps negative Ca/Mg, Al/Mg, Fe/Mg to 0; enforces Si/Mg > 0
 *  - enforces middle-mantle feasibility Al/Mg >= (2/3)Ca/Mg by clamping Ca/Mg down if needed
 *  - enforces global minimum Si/Mg needed to host fixed Ca/Al silicates across all mantle layers:
 *        Si/Mg >= max( 2Ca/Mg + 1.5Al/Mg , 1.5Al/Mg , Ca/Mg )
 *  - attempts solve; if it fails, optionally retries with Ca=Al=0
 *  - returns any adjustments in 'warning'
 */
bool compute_all_mantle_fractions(double CaMg,
                                 double SiMg,
                                 double AlMg,
                                 double FeMg,
                                 vector<double> &upper_out,
                                 vector<double> &middle_out,
                                 vector<double> &lower_out,
                                 string &warning)
{
    stringstream warn;
    warning.clear();

    double R_Ca = CaMg;
    double R_Si = SiMg;
    double R_Al = AlMg;
    double R_Fe = FeMg;

    if (R_Ca < 0.0) { warn << "Ca/Mg < 0; set to 0. "; R_Ca = 0.0; }
    if (R_Al < 0.0) { warn << "Al/Mg < 0; set to 0. "; R_Al = 0.0; }
    if (R_Fe < 0.0) { warn << "Fe/Mg < 0; set to 0. "; R_Fe = 0.0; }
    if (R_Si <= 0.0) { warn << "Si/Mg <= 0 is unphysical; set to 0.1. "; R_Si = 0.1; }

    // Middle mantle feasibility: require R_Al >= (2/3) R_Ca
    if (R_Ca > 0.0 && R_Al < (2.0/3.0)*R_Ca - EPS) {
        double R_Ca_old = R_Ca;
        R_Ca = 1.5 * R_Al;
        warn << "Al/Mg too low to host requested Ca/Mg in middle mantle (need Al/Mg >= 2/3 Ca/Mg). "
             << "Clamping Ca/Mg from " << R_Ca_old << " to " << R_Ca << ". ";
    }

    // Global minimum Si/Mg to host fixed Ca/Al phases in all mantle layers
    const double Si_min_upper  = 2.0*R_Ca;
    const double Si_min_middle = 1.5*R_Al;
    const double Si_min_lower  = R_Ca;

    const double Si_min_global =
        (std::max)(Si_min_upper, (std::max)(Si_min_middle, Si_min_lower));

    if (R_Si < Si_min_global - EPS) {
        double R_Si_old = R_Si;
        R_Si = Si_min_global;
        warn << "Si/Mg too low to host fixed Ca/Al phases in all mantle layers. "
             << "Need Si/Mg >= " << Si_min_global
             << ". Clamping Si/Mg from " << R_Si_old << " to " << R_Si << ". ";
    }

    // First attempt
    bool ok_upper  = compute_upper_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, upper_out);
    bool ok_middle = compute_middle_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, middle_out);
    bool ok_lower  = compute_lower_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, lower_out);

    if (ok_upper && ok_middle && ok_lower) {
        warning = warn.str();
        return true;
    }

    // Optional fallback while middle/lower are still being updated
    warn << "Mantle stoichiometry failed after clamping. "
         << "ok_upper=" << ok_upper << " ok_middle=" << ok_middle << " ok_lower=" << ok_lower << ". "
         << "Retrying with Ca/Mg=0 and Al/Mg=0. ";

    R_Ca = 0.0;
    R_Al = 0.0;

    ok_upper  = compute_upper_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, upper_out);
    ok_middle = compute_middle_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, middle_out);
    ok_lower  = compute_lower_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, lower_out);

    if (ok_upper && ok_middle && ok_lower) {
        warning = warn.str();
        return true;
    }

    warning = warn.str() + " Mantle composition could not be represented even after fallbacks.";
    return false;
}

  // ---- knobs you can tune: mass fraction of Mg-endmember ----
  static const double x_Fo_mass  = 0.5;  // Fo fraction in Fo–Fay mix

  
  // -------------------- LIBRARY OF MIXTURES --------------------
  // -------------------- Fo–Fay 50/50 mixture --------------------
  static vector<EOS*> comps_FoFay{Fo, Ice_Seager};
  static vector<double> x_FoFay{x_Fo_mass, 1.0 - x_Fo_mass};
  DEFINE_IDEAL_MIX_WRAPPERS(FoFay)

  // -------------------- Upper mantle, Olivine Region Mixture --------------------
  static vector<EOS*> comps_OlMix{Forsterite, Fayalite, Enstatite, Ferrosilite, Stishovite, Periclase, Wustite, Diopside, Hedenbergite, Pyrope, Almandine, Spinel, Hercynite};
  static vector<double> x_OlMix(13, 0.0);  // will be filled at runtime
  DEFINE_IDEAL_MIX_WRAPPERS(OlMix)

    // -------------------- Upper mantle, Olivine Region Mixture with Coesite--------------------
  static vector<EOS*> comps_OlMixCoes{Forsterite, Fayalite, Enstatite, Ferrosilite, Coesite, Periclase, Wustite, Diopside, Hedenbergite, Pyrope, Almandine, Spinel, Hercynite};
  static vector<double> x_OlMixCoes(13, 0.0);  // will be filled at runtime
  DEFINE_IDEAL_MIX_WRAPPERS(OlMixCoes)

  // -------------------- Middle mantle, Wadsleyite Region Mixture --------------------
  static vector<EOS*> comps_WdsMix{Mg_Wadsleyite , Fe_Wadsleyite, Stishovite, Periclase, Wustite, Pyrope, Almandine, Grossular, Mg_Majorite, Fe_Akimotoite};
  static vector<double> x_WdsMix(10,0);
  DEFINE_IDEAL_MIX_WRAPPERS(WdsMix)

  // -------------------- Middle mantle, Ringwoodite Region Mixture --------------------
  static vector<EOS*> comps_RwdMix{Mg_Ringwoodite, Fe_Ringwoodite, Stishovite, Periclase, Wustite, Pyrope, Almandine, Grossular, Mg_Majorite, Fe_Akimotoite};
  static vector<double> x_RwdMix(10,0);
  DEFINE_IDEAL_MIX_WRAPPERS(RwdMix)
  
  // -------------------- Lower mantle, Bridgmanite Region Mixture --------------------
  static vector<EOS*> comps_BrgMix{Mg_Bridgmanite, Fe_Bridgmanite, Al_Bridgmanite, Stishovite, Periclase, Wustite, Ca_Perovskite};
  static vector<double> x_BrgMix(7,0);
  DEFINE_IDEAL_MIX_WRAPPERS(BrgMix)

  // -------------------- Lower mantle, Post-Perovskite Region Mixture --------------------
  static vector<EOS*> comps_PPvMix{Mg_Post_Perovskite, Fe_Post_Perovskite, Al_Post_Perovskite, Stishovite, Periclase, Wustite, Ca_Perovskite};
  static vector<double> x_PPvMix(7,0);
  DEFINE_IDEAL_MIX_WRAPPERS(PPvMix)
 
  // -------------------- Rock-Water Mix --------------------
  static vector<EOS*> comps_RockWatMix{Fo, H2O_AQUA};
  static vector<double> x_RockWatMix{0.9269, 0.0731};
  DEFINE_IDEAL_MIX_WRAPPERS(RockWatMix)

  // -------------------- Atm Mix --------------------
  static vector<EOS*> comps_AtmMix{Gas_hhe, watervapor};
  static vector<double> x_AtmMix{0.5, 0.5};
  DEFINE_IDEAL_MIX_WRAPPERS(AtmMix)


  // Sets mantle mixtures from user defined ratios using compute_all_mantle_fractions()
  bool set_mantle_mixtures_from_ratios(double CaMg,
                                     double SiMg,
                                     double AlMg,
                                     double FeMg,
                                     string &warning)
  {
    vector<double> upper_out;
    vector<double> middle_out;
    vector<double> lower_out;
    warning.clear();

    bool ok = compute_all_mantle_fractions(
        CaMg, SiMg, AlMg, FeMg,
        upper_out, middle_out, lower_out,
        warning
    );

    if (!ok) {
        return false;  // caller can decide to abort or fall back
    }

    // Copy into the static vectors used by the mix wrappers
    if (upper_out.size() == x_OlMix.size())
    {
        x_OlMix = upper_out;
        x_OlMixCoes = upper_out;
    }
    else
        return false;

    // index 4 is the SiO2 slot in {Fo,Fa,En,Fs,SiO2,MgO,FeO,Di,Hd,Py,Alm}
    constexpr size_t IDX_SIO2 = 4;
    constexpr double SIO2_ACTIVE_EPS = 1e-10;
    g_upper_mantle_has_free_sio2 = (x_OlMix[IDX_SIO2] > SIO2_ACTIVE_EPS);

    if (middle_out.size() == x_WdsMix.size()){
        x_WdsMix = middle_out;
        x_RwdMix = middle_out;}
    else
        return false;

    if (lower_out.size() == x_BrgMix.size()){
        x_BrgMix = lower_out;
        x_PPvMix = lower_out;}
    else
        return false;

    return true;
  }
} // namespace Mixing
