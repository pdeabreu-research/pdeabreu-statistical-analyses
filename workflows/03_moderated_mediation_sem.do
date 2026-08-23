/****************************************************************************************
Moderated mediation and serial indirect effects using structural equation models
Author: Pedro De Abreu
Public portfolio version: sanitized and generalized from prior research workflows.

Methods demonstrated
  - factorial interaction models and margins
  - parallel multiple-mediator SEM
  - conditional indirect effects via nlcom
  - nonparametric bootstrap inference for conditional indirect effects
  - serial moderated mediation

The script generates synthetic data by default. No participant IDs, collaborator names,
local file paths, unpublished stimuli, or proprietary data are included.
****************************************************************************************/

version 18
clear all
set more off
set seed 20260823

capture mkdir "outputs"

* ------------------------------------------------------------------------------
* 1. Synthetic experimental data
* ------------------------------------------------------------------------------

set obs 900
gen long id = _n
gen double x = rnormal()
gen byte w = runiform() > .50
gen double covariate = rnormal()

summarize x, meanonly
gen double x_c = x - r(mean)
gen double xw = x_c * w

* Two parallel mediators and two downstream outcomes.
gen double m1 = .30*x_c + .18*w + .24*xw + .10*covariate + rnormal()
gen double m2 = -.18*x_c + .15*w - .20*xw + .08*covariate + rnormal()
gen double y = .28*m1 - .22*m2 + .10*x_c + .08*w + .06*xw + .08*covariate + rnormal()
gen double y_serial = .30*y + .16*m1 - .10*m2 + .06*x_c + rnormal()

* ------------------------------------------------------------------------------
* 2. Factorial interaction and marginal effects
* ------------------------------------------------------------------------------

regress y c.x_c##i.w covariate, vce(robust)
margins w, dydx(x_c)
margins, at(x_c=(-2 -1 0 1 2) w=(0 1))
marginsplot, xdimension(x_c) noci name(interaction_plot, replace)
graph export "outputs/interaction_margins.png", replace width(1800)

* ------------------------------------------------------------------------------
* 3. Parallel moderated mediation
* ------------------------------------------------------------------------------
* Moderation operates on the x -> mediator paths. Conditional indirect effects
* are evaluated at w = 0 and w = 1.

sem ///
    (m1 <- x_c w xw covariate) ///
    (m2 <- x_c w xw covariate) ///
    (y  <- m1 m2 x_c w xw covariate), method(mlmv)

estat gof, stats(all)

* Mediator 1 conditional indirect effects
nlcom (ind_m1_w0: (_b[m1:x_c] + 0*_b[m1:xw]) * _b[y:m1]) ///
      (ind_m1_w1: (_b[m1:x_c] + 1*_b[m1:xw]) * _b[y:m1]) ///
      (diff_m1:   ((_b[m1:x_c] + 1*_b[m1:xw]) - ///
                   (_b[m1:x_c] + 0*_b[m1:xw])) * _b[y:m1])

* Mediator 2 conditional indirect effects
nlcom (ind_m2_w0: (_b[m2:x_c] + 0*_b[m2:xw]) * _b[y:m2]) ///
      (ind_m2_w1: (_b[m2:x_c] + 1*_b[m2:xw]) * _b[y:m2]) ///
      (diff_m2:   ((_b[m2:x_c] + 1*_b[m2:xw]) - ///
                   (_b[m2:x_c] + 0*_b[m2:xw])) * _b[y:m2])

capture program drop boot_parallel_modmed
program define boot_parallel_modmed, rclass
    quietly sem ///
        (m1 <- x_c w xw covariate) ///
        (m2 <- x_c w xw covariate) ///
        (y  <- m1 m2 x_c w xw covariate), method(mlmv)

    return scalar m1_w0 = (_b[m1:x_c] + 0*_b[m1:xw]) * _b[y:m1]
    return scalar m1_w1 = (_b[m1:x_c] + 1*_b[m1:xw]) * _b[y:m1]
    return scalar m1_diff = ((_b[m1:x_c] + 1*_b[m1:xw]) - (_b[m1:x_c] + 0*_b[m1:xw])) * _b[y:m1]

    return scalar m2_w0 = (_b[m2:x_c] + 0*_b[m2:xw]) * _b[y:m2]
    return scalar m2_w1 = (_b[m2:x_c] + 1*_b[m2:xw]) * _b[y:m2]
    return scalar m2_diff = ((_b[m2:x_c] + 1*_b[m2:xw]) - (_b[m2:x_c] + 0*_b[m2:xw])) * _b[y:m2]
end

bootstrap ///
    r(m1_w0) r(m1_w1) r(m1_diff) ///
    r(m2_w0) r(m2_w1) r(m2_diff), ///
    reps(5000) seed(20260823) nodots: boot_parallel_modmed
estat bootstrap, all

* ------------------------------------------------------------------------------
* 4. Serial moderated mediation
* ------------------------------------------------------------------------------
* x -> m1 -> y -> y_serial, with moderation on x -> m1.

sem ///
    (m1       <- x_c w xw covariate) ///
    (y        <- m1 m2 x_c w xw covariate) ///
    (y_serial <- y m1 m2 x_c w xw covariate), method(mlmv)

nlcom ///
    (serial_w0: (_b[m1:x_c] + 0*_b[m1:xw]) * _b[y:m1] * _b[y_serial:y]) ///
    (serial_w1: (_b[m1:x_c] + 1*_b[m1:xw]) * _b[y:m1] * _b[y_serial:y]) ///
    (serial_diff: (((_b[m1:x_c] + 1*_b[m1:xw]) - ///
                    (_b[m1:x_c] + 0*_b[m1:xw])) * _b[y:m1] * _b[y_serial:y]))

capture program drop boot_serial_modmed
program define boot_serial_modmed, rclass
    quietly sem ///
        (m1       <- x_c w xw covariate) ///
        (y        <- m1 m2 x_c w xw covariate) ///
        (y_serial <- y m1 m2 x_c w xw covariate), method(mlmv)

    return scalar serial_w0 = (_b[m1:x_c] + 0*_b[m1:xw]) * _b[y:m1] * _b[y_serial:y]
    return scalar serial_w1 = (_b[m1:x_c] + 1*_b[m1:xw]) * _b[y:m1] * _b[y_serial:y]
    return scalar serial_diff = (((_b[m1:x_c] + 1*_b[m1:xw]) - (_b[m1:x_c] + 0*_b[m1:xw])) * _b[y:m1] * _b[y_serial:y])
end

bootstrap r(serial_w0) r(serial_w1) r(serial_diff), ///
    reps(5000) seed(20260824) nodots: boot_serial_modmed
estat bootstrap, all

di as result "Complete. Outputs written to ./outputs/"
