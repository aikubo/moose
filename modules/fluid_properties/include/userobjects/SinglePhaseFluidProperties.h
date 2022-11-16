//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FluidProperties.h"

/**
 * Adds AD versions of each fluid property. These functions use the Real versions of these methods
 * to compute the AD variables complete with derivatives. Typically, these do not need to be
 * overriden in derived classes.
 */
#define propfuncAD(want, prop1, prop2)                                                             \
  virtual DualReal want##_from_##prop1##_##prop2(const DualReal & p1, const DualReal & p2) const   \
  {                                                                                                \
    Real x = 0;                                                                                    \
    Real raw1 = p1.value();                                                                        \
    Real raw2 = p2.value();                                                                        \
    Real dxd1 = 0;                                                                                 \
    Real dxd2 = 0;                                                                                 \
    want##_from_##prop1##_##prop2(raw1, raw2, x, dxd1, dxd2);                                      \
                                                                                                   \
    DualReal result = x;                                                                           \
    result.derivatives() = p1.derivatives() * dxd1 + p2.derivatives() * dxd2;                      \
    return result;                                                                                 \
  }                                                                                                \
                                                                                                   \
  virtual void want##_from_##prop1##_##prop2(const DualReal & prop1,                               \
                                             const DualReal & prop2,                               \
                                             DualReal & val,                                       \
                                             DualReal & d##want##d1,                               \
                                             DualReal & d##want##d2) const                         \
  {                                                                                                \
    fluidPropError(name(), ": ", __PRETTY_FUNCTION__, " derivative derivatives not implemented."); \
    Real dummy, tmp1, tmp2;                                                                        \
    val = want##_from_##prop1##_##prop2(prop1, prop2);                                             \
    want##_from_##prop1##_##prop2(prop1.value(), prop2.value(), dummy, tmp1, tmp2);                \
    d##want##d1 = tmp1;                                                                            \
    d##want##d2 = tmp2;                                                                            \
  }

/**
 * Adds function definitions with not implemented error. These functions should be overriden in
 * derived classes where required. AD versions are constructed automatically using propfuncAD.
 */
#define propfunc(want, prop1, prop2)                                                               \
  virtual Real want##_from_##prop1##_##prop2(Real, Real) const                                     \
  {                                                                                                \
    mooseError(name(), ": ", __PRETTY_FUNCTION__, " not implemented.");                            \
  }                                                                                                \
                                                                                                   \
  virtual void want##_from_##prop1##_##prop2(                                                      \
      Real prop1, Real prop2, Real & val, Real & d##want##d1, Real & d##want##d2) const            \
  {                                                                                                \
    fluidPropError(name(), ": ", __PRETTY_FUNCTION__, " derivatives not implemented.");            \
    d##want##d1 = 0;                                                                               \
    d##want##d2 = 0;                                                                               \
    val = want##_from_##prop1##_##prop2(prop1, prop2);                                             \
  }                                                                                                \
                                                                                                   \
  propfuncAD(want, prop1, prop2)

/**
 * Adds Real declarations of functions that have a default implementation.
 * Important: properties declared using this macro must be defined in SinglePhaseFluidProperties.C.
 * AD versions are constructed automatically using propfuncAD.
 */
#define propfuncWithDefault(want, prop1, prop2)                                                    \
  virtual Real want##_from_##prop1##_##prop2(Real, Real) const;                                    \
  virtual void want##_from_##prop1##_##prop2(                                                      \
      Real prop1, Real prop2, Real & val, Real & d##want##d1, Real & d##want##d2) const;           \
                                                                                                   \
  propfuncAD(want, prop1, prop2)

/**
 * Common class for single phase fluid properties
 */
class SinglePhaseFluidProperties : public FluidProperties
{
public:
  static InputParameters validParams();

  SinglePhaseFluidProperties(const InputParameters & parameters);
  virtual ~SinglePhaseFluidProperties();

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
  // clang-format off

  /**
   * @brief Compute a fluid property given for the state defined by two given properties.
   *
   * For all functions, the first two arguments are the given properties that define the fluid
   * state.  For the two-argument variants, the desired property is the return value.
   * The five-argument variants also provide partial derivatives dx/da and dx/db where x is the
   * desired property being computed, a is the first given property, and b is the second given
   * property.  The desired property, dx/da, and dx/db are stored into the 3rd, 4th, and 5th
   * arguments respectively.
   *
   * Properties/parameters used in these function are listed below with their units:
   *
   * @begincode
   * p      pressure [Pa]
   * T      temperature [K]
   * e      specific internal energy [J/kg]
   * v      specific volume [m^3/kg]
   * rho    density [kg/m^3]
   * h      specific enthalpy [J/kg]
   * s      specific entropy [J/(kg*K)]
   * mu     viscosity [Pa*s]
   * k      thermal conductivity [W/(m*K)]
   * c      speed of sound [m/s]
   * cp     constant-pressure specific heat [J/K]
   * cv     constant-volume specific heat [J/K]
   * beta   volumetric thermal expansion coefficient [1/K]
   * g      Gibbs free energy [J]
   * pp_sat partial pressure at saturation [Pa]
   * gamma  Adiabatic ratio (cp/cv) [-]
   * @endcode
   *
   * As an example:
   *
   * @begincode
   * // calculate pressure given specific vol and energy:
   * auto pressure = your_fluid_properties_object.p_from_v_e(specific_vol, specific_energy);
   *
   * // or use the derivative variant:
   * Real dp_dv = 0; // derivative will be stored into here
   * Real dp_de = 0; // derivative will be stored into here
   * your_fluid_properties_object.p_from_v_e(specific_vol, specific_energy, pressure, dp_dv, dp_de);
   * @endcode
   *
   * Automatic differentiation (AD) support is provided through x_from_a_b(DualReal a, DualReal b) and
   * x_from_a_b(DualReal a, DualReal b, DualReal x, DualReal dx_da, DualReal dx_db) versions of the
   * functions where a and b must be ADReal/DualNumber's calculated using all AD-supporting values:
   *
   * @begincode
   * auto v = 1/rho; // rho must be an AD non-linear variable.
   * auto e = rhoE/rho - vel_energy; // rhoE and vel_energy must be AD variables/numbers also.
   * auto pressure = your_fluid_properties_object.p_from_v_e(v, e);
   * // pressure now contains partial derivatives w.r.t. all degrees of freedom
   * @endcode
   */
  ///@{
  propfunc(p, v, e)
  propfunc(T, v, e)
  propfunc(c, v, e)
  propfunc(cp, v, e)
  propfunc(cv, v, e)
  propfunc(mu, v, e)
  propfunc(k, v, e)
  propfuncWithDefault(s, v, e)
  propfunc(s, h, p)
  propfunc(rho, p, s)
  propfunc(e, v, h)
  propfuncWithDefault(s, p, T)
  propfunc(pp_sat, p, T)
  propfunc(mu, rho, T)
  propfunc(k, rho, T)
  propfuncWithDefault(c, p, T)
  propfuncWithDefault(cp, p, T)
  propfuncWithDefault(cv, p, T)
  propfuncWithDefault(mu, p, T)
  propfuncWithDefault(k, p, T)
  propfunc(rho, p, T)
  propfunc(e, p, rho)
  propfunc(e, T, v)
  propfunc(p, T, v)
  propfunc(h, T, v)
  propfunc(s, T, v)
  propfunc(cv, T, v)
  propfunc(h, p, T)
  propfunc(g, v, e)
  propfuncWithDefault(p, h, s)
  propfunc(T, h, p)  // temporary, until uniformization
  propfuncWithDefault(T, p, h)
  propfuncWithDefault(beta, p, T)
  propfuncWithDefault(v, p, T)
  propfuncWithDefault(e, p, T)
  propfuncWithDefault(gamma, v, e)
  propfuncWithDefault(gamma, p, T)
  ///@}

  // clang-format on

#undef propfunc
#undef propfuncWithDefault
#undef propfuncAD

      /**
       * Fluid name
       * @return string representing fluid name
       */
      virtual std::string fluidName() const;

  /**
   * Molar mass [kg/mol]
   * @return molar mass
   */
  virtual Real molarMass() const;

  /**
   * Critical pressure
   * @return critical pressure (Pa)
   */
  virtual Real criticalPressure() const;

  /**
   * Critical temperature
   * @return critical temperature (K)
   */
  virtual Real criticalTemperature() const;

  /**
   * Critical density
   * @return critical density (kg/m^3)
   */
  virtual Real criticalDensity() const;

  /**
   * Critical specific internal energy
   * @return specific internal energy (J/kg)
   */
  virtual Real criticalInternalEnergy() const;

  /**
   * Triple point pressure
   * @return triple point pressure (Pa)
   */
  virtual Real triplePointPressure() const;

  /**
   * Triple point temperature
   * @return triple point temperature (K)
   */
  virtual Real triplePointTemperature() const;

  /**
   * Specific internal energy from temperature and specific volume
   *
   * @param[in] T     temperature
   * @param[in] v     specific volume
   */
  virtual Real e_spndl_from_v(Real v) const;

  /**
   * Specific internal energy from temperature and specific volume
   *
   * @param[in] T     temperature
   * @param[in] v     specific volume
   */
  virtual void v_e_spndl_from_T(Real T, Real & v, Real & e) const;

  /**
   * Vapor pressure. Used to delineate liquid and gas phases.
   * Valid for temperatures between the triple point temperature
   * and the critical temperature
   *
   * @param T fluid temperature (K)
   * @param[out] saturation pressure (Pa)
   * @param[out] derivative of saturation pressure wrt temperature (Pa/K)
   */
  virtual Real vaporPressure(Real T) const;
  virtual void vaporPressure(Real T, Real & psat, Real & dpsat_dT) const;
  DualReal vaporPressure(const DualReal & T) const;

  /**
   * Vapor temperature. Used to delineate liquid and gas phases.
   * Valid for pressures between the triple point pressure
   * and the critical pressure
   *
   * @param p fluid pressure (Pa)
   * @param[out] saturation temperature (K)
   * @param[out] derivative of saturation temperature wrt pressure
   */
  virtual Real vaporTemperature(Real p) const;
  virtual void vaporTemperature(Real p, Real & Tsat, Real & dTsat_dp) const;
  DualReal vaporTemperature(const DualReal & p) const;

  /**
   * Henry's law coefficients for dissolution in water
   * @return Henry's constant coefficients
   */
  virtual std::vector<Real> henryCoefficients() const;

  virtual void v_e_from_p_T(Real p, Real T, Real & v, Real & e) const;
  virtual void v_e_from_p_T(Real p,
                            Real T,
                            Real & v,
                            Real & dv_dp,
                            Real & dv_dT,
                            Real & e,
                            Real & de_dp,
                            Real & de_dT) const;

  /**
   * Combined methods. These methods are particularly useful for the PorousFlow
   * module, where density and viscosity are typically both computed everywhere.
   * The combined methods allow the most efficient means of calculating both
   * properties, especially where rho(p, T) and mu(rho, T). In this case, an
   * extra density calculation would be required to calculate mu(p, T). All
   * property names are described above.
   */
  virtual void rho_mu_from_p_T(Real p, Real T, Real & rho, Real & mu) const;
  virtual void rho_mu_from_p_T(Real p,
                               Real T,
                               Real & rho,
                               Real & drho_dp,
                               Real & drho_dT,
                               Real & mu,
                               Real & dmu_dp,
                               Real & dmu_dT) const;
  virtual void
  rho_mu_from_p_T(const DualReal & p, const DualReal & T, DualReal & rho, DualReal & mu) const;

  virtual void rho_e_from_p_T(Real p,
                              Real T,
                              Real & rho,
                              Real & drho_dp,
                              Real & drho_dT,
                              Real & e,
                              Real & de_dp,
                              Real & de_dT) const;

  /**
   * Determines (p,T) from (v,e) using Newton Solve in 2D
   * Useful for conversion between different sets of state variables
   *
   * @param[in] v specific volume (m^3 / kg)
   * @param[in] e specific internal energy (J / kg)
   * @param[in] p0 initial guess for pressure (Pa / kg)
   * @param[in] T0 initial guess for temperature (K)
   * @param[out] fluid pressure (Pa / kg)
   * @param[out] Temperature (K)
   */
  virtual void p_T_from_v_e(const Real & v,
                            const Real & e,
                            const Real & p0,
                            const Real & T0,
                            Real & p,
                            Real & T,
                            bool & conversion_succeeded) const;

  /**
   * Determines (p,T) from (v,h) using Newton Solve in 2D
   * Useful for conversion between different sets of state variables
   *
   * @param[in] v specific volume (m^3 / kg)
   * @param[in] h specific enthalpy (J / kg)
   * @param[in] p0 initial guess for pressure (Pa / kg)
   * @param[in] T0 initial guess for temperature (K)
   * @param[out] fluid pressure (Pa / kg)
   * @param[out] Temperature (K)
   */
  virtual void p_T_from_v_h(const Real & v,
                            const Real & h,
                            const Real & p0,
                            const Real & T0,
                            Real & p,
                            Real & T,
                            bool & conversion_succeeded) const;
  /**
   * Determines (p,T) from (h,s) using Newton Solve in 2D
   * Useful for conversion between different sets of state variables
   *
   * @param[in] h specific enthalpy (J / kg)
   * @param[in] s specific entropy (J/K*kg)
   * @param[in] p0 initial guess for pressure (Pa / kg)
   * @param[in] T0 initial guess for temperature (K)
   * @param[out] fluid pressure (Pa / kg)
   * @param[out] Temperature (K)
   */
  virtual void p_T_from_h_s(const Real & h,
                            const Real & s,
                            const Real & p0,
                            const Real & T0,
                            Real & p,
                            Real & T,
                            bool & conversion_succeeded) const;

  template <typename T, typename Functor, typename ADFunctor>
  static void xyDerivatives(const T x,
                            const T & y,
                            T & z,
                            T & dz_dx,
                            T & dz_dy,
                            const Functor & z_from_x_y,
                            const ADFunctor & z_from_x_y_ad);

  /**
   * Newton's method may be used to convert between variable sets
   * _tolerance, _T_initial_guess, and _p_initial_guess are the parameters for these
   * iterative solves
   */
  const Real _tolerance;
  const Real _T_initial_guess;
  const Real _p_initial_guess;

private:
  template <typename... Args>
  void fluidPropError(Args... args) const
  {
    if (_allow_imperfect_jacobians)
      mooseDoOnce(mooseWarning(std::forward<Args>(args)...));
    else
      mooseError(std::forward<Args>(args)...);
  }

  template <typename T>
  static void assignQoIs(const T & /*x*/,
                         const T & /*y*/,
                         T & /*z*/,
                         T & /*dz_dx*/,
                         T & /*dz_dy*/,
                         const FPDualReal & /*raw_z*/);

  template <typename T>
  static void checkDerivatives(const T & /*a*/, const T & /*b*/);
};

#pragma GCC diagnostic pop

template <typename T>
inline void
SinglePhaseFluidProperties::checkDerivatives(const T & a, const T & b)
{
  mooseAssert(MooseUtils::relativeFuzzyEqual(a.value(), b.value(), 1e-2), "Must be equal");
  mooseAssert(a.derivatives().size() == b.derivatives().size(), "Must be equal");
  const auto & a_derivs = a.derivatives().nude_data();
  const auto & b_derivs = b.derivatives().nude_data();
  for (const auto i : make_range(a.derivatives().size()))
    mooseAssert(MooseUtils::relativeFuzzyEqual(a_derivs[i], b_derivs[i], 1e-2), "Must be equal");
}

template <>
inline void
SinglePhaseFluidProperties::checkDerivatives(const Real &, const Real &)
{
}

template <typename T>
void
SinglePhaseFluidProperties::assignQoIs(
    const T & x, const T & y, T & z, T & dz_dx, T & dz_dy, const FPDualReal & raw_z)
{
  z = raw_z.value();
  dz_dx = raw_z.derivatives()[0];
  dz_dy = raw_z.derivatives()[1];

  z.derivatives() = dz_dx.value() * x.derivatives() + dz_dy.value() * y.derivatives();
  dz_dx.derivatives() = x.derivatives();
  dz_dy.derivatives() = y.derivatives();
}

template <>
inline void
SinglePhaseFluidProperties::assignQoIs(const Real & /*x*/,
                                       const Real & /*y*/,
                                       Real & z,
                                       Real & dz_dx,
                                       Real & dz_dy,
                                       const FPDualReal & raw_z)
{
  z = raw_z.value();
  dz_dx = raw_z.derivatives()[0];
  dz_dy = raw_z.derivatives()[1];
}

template <typename T, typename Functor, typename ADFunctor>
void
SinglePhaseFluidProperties::xyDerivatives(const T x,
                                          const T & y,
                                          T & z,
                                          T & dz_dx,
                                          T & dz_dy,
                                          const Functor &
#ifdef DEBUG
                                              z_from_x_y
#endif
                                          ,
                                          const ADFunctor & z_from_x_y_ad)
{
  FPDualReal raw_x = MetaPhysicL::raw_value(x);
  Moose::derivInsert(raw_x.derivatives(), 0, 1);
  FPDualReal raw_y = MetaPhysicL::raw_value(y);
  Moose::derivInsert(raw_y.derivatives(), 1, 1);
  const FPDualReal raw_z = z_from_x_y_ad(raw_x, raw_y);
  assignQoIs(x, y, z, dz_dx, dz_dy, raw_z);

#ifdef DEBUG
  const auto straight_z = z_from_x_y(x, y);
  checkDerivatives(z, straight_z);
  static constexpr Real perturbation_factor = 1 + 1e-8;
  const auto perturbed_z_x = z_from_x_y(perturbation_factor * x, y);
  const auto perturbed_z_y = z_from_x_y(x, perturbation_factor * y);
  const auto Jx = (perturbed_z_x - z) / (1e-8 * x);
  const auto Jy = (perturbed_z_y - z) / (1e-8 * y);
  if (!MooseUtils::relativeFuzzyEqual(Jx, dz_dx, 1e-2))
  {
    const FPDualReal raw_z = z_from_x_y_ad(raw_x, raw_y);
    assignQoIs(x, y, z, dz_dx, dz_dy, raw_z);
    ::mooseError("Bad x Jacobian");
  }
  if (!MooseUtils::relativeFuzzyEqual(Jy, dz_dy, 1e-2))
  {
    const FPDualReal raw_z = z_from_x_y_ad(raw_x, raw_y);
    assignQoIs(x, y, z, dz_dx, dz_dy, raw_z);
    ::mooseError("Bad y Jacobian");
  }
#endif
}
