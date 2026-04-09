// mixing.h
#ifndef MIXING_H
#define MIXING_H

#include "EOS.h"

// Keep mixing tools in a namespace so they don't collide
namespace Mixing {

  // ---------- General helpers (any number of components) ----------

  // Ideal mixture density:
  // - P_cgs in microbar (same as EOS::density)
  // - T in K
  // - components: EOS* of each pure phase
  // - x_mass: mass fractions of each component
  double density_ideal_mixture(double P_cgs, double T,
                               const vector<EOS*> &components,
                               const vector<double> &x_mass, double rho_guess = 1.0);

  // Ideal mixture adiabatic gradient dT/dP|S:
  // - P in cgs (same as EOS::dTdP_S)
  // - returns K/GPa
  double dTdP_S_ideal_mixture(double P_cgs, double T,
                              const vector<EOS*> &components,
                              const vector<double> &x_mass, double rho_guess = 1.0);

  // Ideal mixture adiabatic gradient dT/dP|S (full rigorous version):
  // - More rigorous formulation accounting for volume expansion
  // - P in cgs (same as EOS::dTdP_S)
  // - returns K/GPa
  double dTdP_S_ideal_mixture_full(double P_cgs, double T,
                                    const vector<EOS*> &components,
                                    const vector<double> &x_mass, double rho_guess = 1.0);

  // ---------- Wrappers for Mixtures ----------
  // one macro to declare the triple definitions for any NAME
  // Density(P,T,rho_guess) in cgs 
  // dT/dP (K/microbar) for Fo+Fay 50/50, for use in EOS::dTdm
  // Mixture dT/dP|S in K/GPa
  #define DECLARE_IDEAL_MIX(NAME)                                      \
    double density_##NAME(double P_cgs, double T, double rho_guess);   \
    double dTdP_##NAME(double P_cgs, double T, double &rho_guess);     \
    double dTdP_S_##NAME(double P_cgs, double T, double &rho_guess);

  //Declaration of Mixtures in Library
  DECLARE_IDEAL_MIX(FoFay)
  DECLARE_IDEAL_MIX(OlMix)
  DECLARE_IDEAL_MIX(OlMixCoes)
  DECLARE_IDEAL_MIX(WdsMix)
  DECLARE_IDEAL_MIX(RwdMix)
  DECLARE_IDEAL_MIX(BrgMix)
  DECLARE_IDEAL_MIX(PPvMix)
  DECLARE_IDEAL_MIX(RockWatMix)
  DECLARE_IDEAL_MIX(AtmMix)

  #undef DECLARE_IDEAL_MIX

  // Convert mantle Mg# (= Mg/(Mg+Fe), molar) to bulk mantle FeO mass fraction
  bool mantle_wFeO_from_MgNumber(
    double CaMg, double SiMg, double AlMg,
    double Mg_number,
    double &mantle_wFeO_out,
    std::string &note_or_error);

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
    std::string &note_or_error);

  // Upper mantle (olivine region) mass fractions
  bool compute_upper_mantle_fractions(double CaMg,
                                      double SiMg,
                                      double AlMg,
                                      double FeMg,
                                      vector<double> &upper_out);

  // Middle mantle (wadsleyite / ringwoodite region) mass fractions
  bool compute_middle_mantle_fractions(double CaMg,
                                      double SiMg,
                                      double AlMg,
                                      double FeMg,
                                      vector<double> &middle_out);

  // Lower mantle (perovskite / post-perovskite region) mass fractions
  bool compute_lower_mantle_fractions(double CaMg,
                                      double SiMg,
                                      double AlMg,
                                      double FeMg,
                                      vector<double> &lower_out);

  // Convenience wrapper: computes all three layers,
  // applies clamping/fallbacks, and returns a warning string
  // describing any adjustments.
  bool compute_all_mantle_fractions(double CaMg,
                                    double SiMg,
                                    double AlMg,
                                    double FeMg,
                                    vector<double> &upper_out,
                                    vector<double> &middle_out,
                                    vector<double> &lower_out,
                                    string &warning);

  // Sets mantle mixtures from user defined ratios using compute_all_mantle_fractions()
  bool set_mantle_mixtures_from_ratios(double CaMg,
                                         double SiMg,
                                         double AlMg,
                                         double FeMg,
                                         string &warning);

  // True when upper mantle mixture has non-negligible free SiO2 (mass fraction > ~1e-10)
  bool upper_mantle_has_free_sio2();

} // namespace Mixing


#endif // MIXING_H
