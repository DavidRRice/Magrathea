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


  // Simple normalization helper
  static void normalize(vector<double> &x)
  {
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
      sum += x[i];

    if (sum <= 0.0) return;

    for (size_t i = 0; i < x.size(); ++i)
      x[i] /= sum;
  }

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
    normalize(x);

    // sanitize guess (avoid absurd values propagating into component solvers)
    double guess = rho_guess;

    double M_mix = 0.0;   // g/mol
    double V_mix = 0.0;   // cm^3/mol

    for (size_t i = 0; i < n; ++i)
    {
      EOS *phase = components[i];
      const double Mi = phase->getmmol();        // g/mol
      if (!gsl_finite(Mi) || Mi <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      const double rho_i = phase->density(P_cgs, T, guess); // g/cm^3
      if (!gsl_finite(rho_i) || rho_i <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      const double v_i = Mi / rho_i;             // cm^3/mol

      M_mix += x[i] * Mi;
      V_mix += x[i] * v_i;
    }

    if (V_mix <= 0.0)
      return numeric_limits<double>::quiet_NaN();

    return M_mix / V_mix; // g/cm^3
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
  double dTdP_S_ideal_mixture(double P_GPa, double T,
                              const vector<EOS*> &components,
                              const vector<double> &x_in, double rho_guess)
  {
    const size_t n = components.size();
    if (n == 0 || x_in.size() != n)
      return numeric_limits<double>::quiet_NaN();

    vector<double> x = x_in;
    normalize(x);

    const double P_cgs = P_GPa * 1.0e10; // microbar

    double inv_grad_sum = 0.0;  // Σ x_i / grad_i using original x_i
    double x_used_sum   = 0.0;  // Σ x_i over phases with a usable gradient
    int    n_thermal    = 0;    // number of phases contributing to gradient

    double guess = rho_guess;

    for (size_t i = 0; i < n; ++i)
    {
      EOS *phase = components[i];

      // Skip EOS that have no thermal parameters: thermal_type == 0
      const int thermal_type = phase->getthermal();
      if (thermal_type == 0)
        continue;

      // We need a reasonable rho_guess for dTdP_S. Use pure-phase density.
      double rho_i = phase->density(P_cgs, T, guess);
      if (!gsl_finite(rho_i) || rho_i <= 0.0)
        return numeric_limits<double>::quiet_NaN();

      const double grad_i = phase->dTdP_S(P_GPa, T, rho_i); // K/GPa

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


  // --- Generic wrapper generator for a given mixture NAME ---
  //
  // Expects:
  //   vector<EOS*> comps_NAME;
  //   vector<double> x_NAME;
  //
  // and generates:
  //   density_NAME(P_cgs, T, rho_guess)
  //   dTdP_S_NAME(P_GPa, T, rho_guess)
  //   dTdP_NAME(P_cgs, T, rho_guess)
  //
  #define DEFINE_IDEAL_MIX_WRAPPERS(NAME)                                       \
  double density_##NAME(double P_cgs, double T, double rho_guess)         \
  {                                                                           \
    return density_ideal_mixture(P_cgs, T, comps_##NAME, x_##NAME, rho_guess);          \
  }                                                                           \
                                                                              \
  double dTdP_S_##NAME(double P_GPa, double T, double &rho_guess)         \
  {                                                                           \
    return dTdP_S_ideal_mixture(P_GPa, T, comps_##NAME, x_##NAME);           \
  }                                                                           \
                                                                              \
  double dTdP_##NAME(double P_cgs, double T, double &rho_guess)               \
  {                                                                           \
    const double P_GPa   = P_cgs / 1.0e10;                                    \
    const double grad_GPa = dTdP_S_ideal_mixture(P_GPa, T,                    \
                                                 comps_##NAME, x_##NAME);     \
    if (!gsl_finite(grad_GPa))                                                \
      return numeric_limits<double>::quiet_NaN();                        \
    return grad_GPa / 1.0e10;                                                 \
  }

  bool mantle_wFeO_from_MgNumber(
  double CaMg, double SiMg, double AlMg,
  double Mg_number,
  double &mantle_wFeO_out,
  std::string &note_or_error)
  {
  note_or_error.clear();
  std::stringstream msg;

  if (!gsl_finite(Mg_number) || Mg_number <= 0.0 || Mg_number >= 1.0) {
    note_or_error = "mantle_Mg_number must be in (0,1).";
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

  const double w_nonFe = wt_fract_S_core + wt_fract_Ni_core;
  if (w_nonFe >= 1.0 - EPS) {
    msg << "Core non-Fe mass fraction (S+Ni)=" << w_nonFe << " must be < 1.";
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
      << ", core non-Fe (S+Ni)=" << w_nonFe << ".";
  note_or_error = msg.str();
  return true;
  }


  
  //---------- MANTLE MIXTURE GENERATOR --------------
  /*
  * upper_out (size 8):
  *   0: Fo          (forsterite, Mg2SiO4)
  *   1: Fayalite    (Fe2SiO4)
  *   2: Enstatite   (Mg2Si2O6)
  *   3: Ferrosilite (Fe2Si2O6)
  *   4: Diopside    (CaMgSi2O6)
  *   5: Hedenbergite( CaFeSi2O6)
  *   6: Pyrope      (Mg3Al2Si3O12)
  *   7: Almandine   (Fe3Al2Si3O12)
  */
bool compute_upper_mantle_fractions(double CaMg,
                                    double SiMg,
                                    double AlMg,
                                    double FeMg,
                                    vector<double> &upper_out)
{
    upper_out.assign(8, 0.0);

    if (CaMg < 0.0 || AlMg < 0.0 || FeMg < 0.0 || SiMg <= 0.0)
        return false;

    const double R_Ca = CaMg;
    const double R_Si = SiMg;
    const double R_Al = AlMg;
    const double R_Fe = FeMg;

    const double X_Mg = 1.0 / (1.0 + R_Fe);
    const double X_Fe = 1.0 - X_Mg;

    // Mineral formula-unit moles (per 1 mole of Mg)
    double n_Cpx = R_Ca;
    double n_Gt  = 0.5 * R_Al;
    double n_Ol  = 1.0 / X_Mg + R_Ca - R_Si; // = 1 + R_Fe + R_Ca - R_Si
    double n_Opx = R_Si - 1.5 * R_Ca - 0.75 * R_Al - 1.0 / (2.0 * X_Mg);  // = R_Si - 1.5*R_Ca -0.75*R_Al -0.5*(1+R_Fe)

    if (!enforce_non_negative(n_Cpx)) return false;
    if (!enforce_non_negative(n_Gt))  return false;
    if (!enforce_non_negative(n_Ol))  return false;
    if (!enforce_non_negative(n_Opx)) return false;

    // Endmember moles
    double n_Fo  = X_Mg * n_Ol;
    double n_Fa  = X_Fe * n_Ol;
    double n_En  = X_Mg * n_Opx;
    double n_Fs  = X_Fe * n_Opx;
    double n_Di  = X_Mg * n_Cpx;
    double n_Hd  = X_Fe * n_Cpx;
    double n_Py  = X_Mg * n_Gt;
    double n_Alm = X_Fe * n_Gt;

    double Ntot = n_Fo + n_Fa + n_En + n_Fs + n_Di + n_Hd + n_Py + n_Alm;

    if (Ntot <= 0.0) return false;

    upper_out[0] = n_Fo  / Ntot;
    upper_out[1] = n_Fa  / Ntot;
    upper_out[2] = n_En  / Ntot;
    upper_out[3] = n_Fs  / Ntot;
    upper_out[4] = n_Di  / Ntot;
    upper_out[5] = n_Hd  / Ntot;
    upper_out[6] = n_Py  / Ntot;
    upper_out[7] = n_Alm / Ntot;

    return true;
}

/*
 * middle_out (size 6):
 *   0: Wadsleyite (Mg) (or Ringwoodite)
 *   1: Fe-Wadsleyite  (or Fe-Ringwoodite)
 *   2: Pyrope
 *   3: Almandine
 *   4: Grossular
 *   5: Mg-Majorite
 */
bool compute_middle_mantle_fractions(double CaMg,
                                     double SiMg,
                                     double AlMg,
                                     double FeMg,
                                     vector<double> &middle_out)
{
    middle_out.assign(6, 0.0);

    if (CaMg < 0.0 || AlMg < 0.0 || FeMg < 0.0 || SiMg <= 0.0)
        return false;

    const double R_Ca = CaMg;
    const double R_Si = SiMg;
    const double R_Al = AlMg;
    const double R_Fe = FeMg;

    double n_Grs = R_Ca / 3.0;
    double n_G   = 0.5 * R_Al - R_Ca / 3.0;
    double n_W   = 1.0 + R_Ca + R_Fe - R_Si;
    double n_Mj  = -0.375 * R_Al - 0.25  * R_Ca - 0.25  * R_Fe + 0.5   * R_Si - 0.25;

    if (!enforce_non_negative(n_Grs)) return false;
    if (!enforce_non_negative(n_G))   return false;
    if (!enforce_non_negative(n_W))   return false;
    if (!enforce_non_negative(n_Mj))  return false;

    // Common Mg# in Wadsleyite + Al-garnet:
    double num =  3.0*R_Al + 2.0*R_Ca + 2.0*R_Fe - 4.0*R_Si + 4.0;
    double den =  3.0*R_Al + 2.0*R_Ca + 4.0*R_Fe - 4.0*R_Si + 4.0;
    if (fabs(den) < EPS) return false;

    double X_Mg = num / den;
    double X_Fe = 1.0 - X_Mg;

    if (X_Mg < -EPS || X_Mg > 1.0 + EPS) return false;
    if (X_Mg < 0.0) X_Mg = 0.0;
    if (X_Mg > 1.0) X_Mg = 1.0;

    double n_Wds    = X_Mg * n_W;
    double n_FeWds  = X_Fe * n_W;
    double n_Py     = X_Mg * n_G;
    double n_Alm    = X_Fe * n_G;
    double n_Gross  = n_Grs;
    double n_MgMaj  = n_Mj;

    double Ntot = n_Wds + n_FeWds + n_Py + n_Alm + n_Gross + n_MgMaj;
    if (Ntot <= 0.0) return false;

    middle_out[0] = n_Wds   / Ntot;
    middle_out[1] = n_FeWds / Ntot;
    middle_out[2] = n_Py    / Ntot;
    middle_out[3] = n_Alm   / Ntot;
    middle_out[4] = n_Gross / Ntot;
    middle_out[5] = n_MgMaj / Ntot;

    return true;
}

/*
 * lower_out (size 6):
 *   0: Mg-Perovskite  (MgSiO3) (or Mg-PPv)
 *   1: Fe-Perovskite  (FeSiO3) (or Fe-PPv)
 *   2: Periclase      (MgO)
 *   3: Wustite        (FeO)
 *   4: Ca-Perovskite  (CaSiO3)
 *   5: Al-Perovskite  (AlAlO3)
 */
bool compute_lower_mantle_fractions(double CaMg,
                                    double SiMg,
                                    double AlMg,
                                    double FeMg,
                                    vector<double> &lower_out)
{
    lower_out.assign(6, 0.0);

    if (CaMg < 0.0 || AlMg < 0.0 || FeMg < 0.0 || SiMg <= 0.0)
        return false;

    const double R_Ca = CaMg;
    const double R_Si = SiMg;
    const double R_Al = AlMg;
    const double R_Fe = FeMg;

    const double X_Mg = 1.0 / (1.0 + R_Fe);
    const double X_Fe = 1.0 - X_Mg;

    double n_CaPv = R_Ca;
    double n_AlPv = 0.5 * R_Al;
    double n_Pv   = R_Si - R_Ca;
    double n_Fp   = 1.0 / X_Mg - n_Pv; // = 1 + R_Fe - (R_Si - R_Ca)

    if (!enforce_non_negative(n_CaPv)) return false;
    if (!enforce_non_negative(n_AlPv)) return false;
    if (!enforce_non_negative(n_Pv))   return false;
    if (!enforce_non_negative(n_Fp))   return false;

    double n_PvMg = X_Mg * n_Pv;
    double n_PvFe = X_Fe * n_Pv;
    double n_FpMg = X_Mg * n_Fp;
    double n_FpFe = X_Fe * n_Fp;

    double Ntot = n_PvMg + n_PvFe + n_FpMg + n_FpFe + n_CaPv + n_AlPv;
    if (Ntot <= 0.0) return false;

    lower_out[0] = n_PvMg / Ntot;
    lower_out[1] = n_PvFe / Ntot;
    lower_out[2] = n_FpMg / Ntot;
    lower_out[3] = n_FpFe / Ntot;
    lower_out[4] = n_CaPv / Ntot;
    lower_out[5] = n_AlPv / Ntot;

    return true;
}

/*
 * Convenience wrapper with fallbacks:
 *  - clamps Si/Mg into its allowed [Si_min, Si_max] range
 *  - if Al/Mg is too small to host Ca/Mg, sets Ca=Al=0
 *  - if something still fails, tries Ca=Al=0 as a last resort
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

    // Copy inputs so we can modify them locally
    double R_Ca = CaMg;
    double R_Si = SiMg;
    double R_Al = AlMg;
    double R_Fe = FeMg;

    // Basic clamping for negative inputs
    if (R_Ca < 0.0) { warn << "Ca/Mg < 0; set to 0. "; R_Ca = 0.0; }
    if (R_Al < 0.0) { warn << "Al/Mg < 0; set to 0. "; R_Al = 0.0; }
    if (R_Fe < 0.0) { warn << "Fe/Mg < 0; set to 0. "; R_Fe = 0.0; }
    if (R_Si <= 0.0) {
        warn << "Si/Mg <= 0 is unphysical; set to 0.1. ";
        R_Si = 0.1;
    }

    // If not enough Al to host Ca, zero both (your requested fallback)
    if (R_Al < (2.0/3.0)*R_Ca && R_Ca > 0.0) {
        warn << "Al/Mg is too low to host Ca/Mg in our mantle model; "
             << "setting Ca/Mg and Al/Mg to 0. ";
        R_Ca = 0.0;
        R_Al = 0.0;
    }

    // Compute allowed Si range for current Ca, Al, Fe:
    //   Si_min from Opx >= 0,
    //   Si_max from Ol/Wds/Fp >= 0.
    auto clamp_Si = [&](double &R_Si_local, double R_Ca_local, double R_Al_local, double R_Fe_local)
    {
        double Si_min = 0.5 * (1.0 + R_Fe_local) + 1.5 * R_Ca_local + 0.75 * R_Al_local;
        double Si_max = 1.0 + R_Fe_local + R_Ca_local;

        if (Si_min > Si_max) {
            // Incompatible; caller will decide on further fallback.
            return std::make_pair(false, std::pair<double,double>(Si_min, Si_max));
        }

        if (R_Si_local < Si_min) {
            warn << "Si/Mg too low; clamped from " << R_Si_local
                 << " to " << Si_min << ". ";
            R_Si_local = Si_min;
        } else if (R_Si_local > Si_max) {
            warn << "Si/Mg too high; clamped from " << R_Si_local
                 << " to " << Si_max << ". ";
            R_Si_local = Si_max;
        }

        return std::make_pair(true, std::pair<double,double>(Si_min, Si_max));
    };

    // First attempt: with (possibly adjusted) Ca, Al, Fe, Si
    {
        auto ok_range = clamp_Si(R_Si, R_Ca, R_Al, R_Fe);
        if (!ok_range.first) {
            // If Si range itself is inconsistent, drop to Ca=Al=0 fallback
            warn << "Input ratios give incompatible Si/Mg bounds; "
                 << "resetting Ca/Mg and Al/Mg to 0. ";
            R_Ca = 0.0;
            R_Al = 0.0;
            // recompute Si bounds for Ca=Al=0
            clamp_Si(R_Si, R_Ca, R_Al, R_Fe);
        }

        bool ok1 = compute_upper_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, upper_out);
        bool ok2 = compute_middle_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, middle_out);
        bool ok3 = compute_lower_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, lower_out);

        if (ok1 && ok2 && ok3) {
            warning = warn.str();
            return true;
        }
    }

    // Last-resort fallback: force Ca=Al=0, clamp Si for that case, and try again
    warn << "Mantle composition still invalid; trying Ca/Mg=0 and Al/Mg=0. ";
    R_Ca = 0.0;
    R_Al = 0.0;
    clamp_Si(R_Si, R_Ca, R_Al, R_Fe);

    bool ok1 = compute_upper_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, upper_out);
    bool ok2 = compute_middle_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, middle_out);
    bool ok3 = compute_lower_mantle_fractions(R_Ca, R_Si, R_Al, R_Fe, lower_out);

    if (ok1 && ok2 && ok3) {
        warning = warn.str();
        return true;
    }

    warning = warn.str() + " Mantle composition could not be represented even after fallbacks.";
    return false;
  }


  // ---- knobs you can tune: molar fraction of Mg-endmember ----
  static const double x_Fo_mol  = 0.5;  // Fo fraction in Fo–Fay mix
  static const double x_Wds_mol = 0.5;  // Mg-Wds fraction in Wds–FeWds mix
  static const double x_Rwd_mol = 0.5;  // Mg-Rwd fraction in Rwd–FeRwd mix

  
  // -------------------- LIBRARY OF MIXTURES --------------------
  // -------------------- Fo–Fay 50/50 mixture --------------------
  static vector<EOS*> comps_FoFay{Fo, Ice_Seager};
  static vector<double> x_FoFay{x_Fo_mol, 1.0 - x_Fo_mol};
  DEFINE_IDEAL_MIX_WRAPPERS(FoFay)

  // -------------------- Upper mantle, Olivine Region Mixture --------------------
  static vector<EOS*> comps_OlMix{Fo, Fayalite, Enstatite, Ferrosilite, Diopside, Hedenbergite, Pyrope, Almandine};
  static vector<double> x_OlMix(8, 0.0);  // will be filled at runtime
  DEFINE_IDEAL_MIX_WRAPPERS(OlMix)

  // -------------------- Upper mantle, Wadsleyite Region Mixture --------------------
  static vector<EOS*> comps_WdsMix{Wds, Fe_Wadsleyite, Pyrope, Almandine, Grossular, Mg_Majorite};
  static vector<double> x_WdsMix(6,0);
  DEFINE_IDEAL_MIX_WRAPPERS(WdsMix)

  // -------------------- Upper mantle, Ringwoodite Region Mixture --------------------
  static vector<EOS*> comps_RwdMix{Rwd, Fe_Ringwoodite, Pyrope, Almandine, Grossular, Mg_Majorite};
  static vector<double> x_RwdMix(6,0);
  DEFINE_IDEAL_MIX_WRAPPERS(RwdMix)
  
  // -------------------- Lower mantle, Bridgmanite Region Mixture --------------------
  static vector<EOS*> comps_BrgMix{Si_Pv, Fe_Perovskite, MgO, B1FeO, Ca_Perovskite, Al_Perovskite};
  static vector<double> x_BrgMix(6,0);
  DEFINE_IDEAL_MIX_WRAPPERS(BrgMix)

  // -------------------- Lower mantle, Post-Perovskite Region Mixture --------------------
  static vector<EOS*> comps_PPvMix{Si_PPv_Sakai, Fe_Post_Perovskite, MgO, B1FeO, Ca_Perovskite, Al_Post_Perovskite};
  static vector<double> x_PPvMix(6,0);
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
        x_OlMix = upper_out;
    else
        return false;

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