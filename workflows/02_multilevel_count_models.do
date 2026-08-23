/****************************************************************************************
Multilevel count models for repeated social-network / engagement data
Author: Pedro De Abreu
Public portfolio version: sanitized and generalized from prior research workflows.

Methods demonstrated
  - repeated-observation data engineering
  - network-normalized predictors and temporal lag construction
  - mixed-effects Poisson count models with account-level random intercepts
  - interaction probing with marginal effects
  - negative-binomial sensitivity analysis

The script generates synthetic data by default. No user handles, platform identifiers,
participant records, or proprietary datasets are included.
****************************************************************************************/

version 18
clear all
set more off
set seed 20260823

capture mkdir "outputs"

* ------------------------------------------------------------------------------
* 1. Synthetic repeated-observation data
* ------------------------------------------------------------------------------
* Replace this section with an authorized import for a real analysis.

set obs 250
gen long account_id = _n
gen double account_re = rnormal(0, .35)
gen double trait_a = rnormal()
gen double trait_b = rnormal()
gen byte group01 = runiform() > .50
gen double age_z = rnormal()
gen double follower_count = round(exp(rnormal(7.5, .70)))
gen double friend_count = round(follower_count * runiform(.10, .80))
gen double network_indegree = abs(rnormal(50, 20))

expand 24
bysort account_id: gen int post_index = _n
bysort account_id: replace account_re = account_re[1]
bysort account_id: replace trait_a = trait_a[1]
bysort account_id: replace trait_b = trait_b[1]
bysort account_id: replace group01 = group01[1]
bysort account_id: replace age_z = age_z[1]
bysort account_id: replace follower_count = follower_count[1]
bysort account_id: replace friend_count = friend_count[1]
bysort account_id: replace network_indegree = network_indegree[1]

gen double positive_affect = rnormal()
gen double negative_affect = rnormal()
gen double time_z = (post_index - 12.5) / 6.9

gen double friend_follower_ratio = friend_count / follower_count
gen double log_network_size = ln(follower_count)
gen double indegree_per_follower = network_indegree / follower_count

* Latent expected engagement rate; kept modest to avoid extreme synthetic counts.
gen double xb = -1.20 + .18*trait_a - .10*trait_b + .16*group01 + ///
    .10*positive_affect - .08*negative_affect + .10*time_z + account_re

gen double mu = exp(xb)
gen int engagement_count = rpoisson(mu)

xtset account_id post_index
gen double lag_engagement = L.engagement_count
replace lag_engagement = 0 if missing(lag_engagement)

gen double rolling_engagement = (L1.engagement_count + L2.engagement_count + L3.engagement_count) / 3
replace rolling_engagement = lag_engagement if missing(rolling_engagement)
replace rolling_engagement = 0 if missing(rolling_engagement)

gen double trait_a_x_group = trait_a * group01
gen double trait_b_x_group = trait_b * group01

* ------------------------------------------------------------------------------
* 2. Descriptives and model sequence
* ------------------------------------------------------------------------------

summarize engagement_count trait_a trait_b friend_follower_ratio ///
    log_network_size indegree_per_follower positive_affect negative_affect

estimates clear

mepoisson engagement_count ///
    friend_follower_ratio log_network_size positive_affect negative_affect ///
    rolling_engagement time_z || account_id:, vce(oim)
estimates store M0

mepoisson engagement_count ///
    trait_a trait_b group01 age_z ///
    friend_follower_ratio log_network_size indegree_per_follower ///
    positive_affect negative_affect rolling_engagement time_z ///
    || account_id:, vce(oim)
estimates store M1

mepoisson engagement_count ///
    c.trait_a##i.group01 c.trait_b##i.group01 age_z ///
    friend_follower_ratio log_network_size indegree_per_follower ///
    positive_affect negative_affect rolling_engagement time_z ///
    || account_id:, vce(oim)
estimates store M2

estimates table M0 M1 M2, b(%9.3f) se(%9.3f) stats(N ll)
estimates stats M0 M1 M2

* ------------------------------------------------------------------------------
* 3. Interaction probing and marginal effects
* ------------------------------------------------------------------------------

margins group01, dydx(trait_a)
margins, at(trait_a=(-2 -1 0 1 2) group01=(0 1)) predict(mu) vsquish
marginsplot, xdimension(trait_a) noci name(traitA_by_group, replace)

graph export "outputs/multilevel_traitA_interaction.png", replace width(1800)

margins group01, dydx(trait_b)
margins, at(trait_b=(-2 -1 0 1 2) group01=(0 1)) predict(mu) vsquish

* ------------------------------------------------------------------------------
* 4. Overdispersion sensitivity check
* ------------------------------------------------------------------------------
* Mixed-effects negative-binomial model provides a useful sensitivity comparison
* when the conditional variance materially exceeds the conditional mean.

menbreg engagement_count ///
    c.trait_a##i.group01 c.trait_b##i.group01 age_z ///
    friend_follower_ratio log_network_size indegree_per_follower ///
    positive_affect negative_affect rolling_engagement time_z ///
    || account_id:
estimates store NB1

estimates table M2 NB1, b(%9.3f) se(%9.3f) stats(N ll)
estimates stats M2 NB1

* ------------------------------------------------------------------------------
* 5. Export a compact coefficient table using built-in collect
* ------------------------------------------------------------------------------

collect clear
collect _r_b _r_se _r_p: mepoisson engagement_count ///
    c.trait_a##i.group01 c.trait_b##i.group01 age_z ///
    friend_follower_ratio log_network_size indegree_per_follower ///
    positive_affect negative_affect rolling_engagement time_z ///
    || account_id:, vce(oim)
collect layout (colname) (result[_r_b _r_se _r_p])
collect export "outputs/multilevel_poisson_results.html", replace

di as result "Complete. Outputs written to ./outputs/"
