options = Options()

options.set('minimizer', 'strategy', 'newton_vanilla')

doUpVariation=False

benchmark='S4Inv1000'

# for model building:
def get_model():
    # Read in and build the model automatically from the histograms in the root file. 
    # This model will contain all shape uncertainties given according to the templates
    # which also includes rate changes according to the alternate shapes.
    # For more info about this model and naming conventuion, see documentation
    # of build_model_from_rootfile.
    
    model = build_model_from_rootfile('inputTheta_WTTinterCRSR.root',include_mc_uncertainties=True)
    #model = build_model_from_rootfile('inputTheta_Wmerged_WTTSR.root',include_mc_uncertainties=True)
    
    # If the prediction histogram is zero, but data is non-zero, teh negative log-likelihood
    # is infinity which causes problems for some methods. Therefore, we set all histogram
    # bin entries to a small, but positive value:
    model.fill_histogram_zerobins(0.0000000001)

    # define what the signal processes are. All other processes are assumed to make up the 
    # 'background-only' model.
    model.set_signal_processes('S4*')



    # Add some lognormal rate uncertainties. The first parameter is the name of the
    # uncertainty (which will also be the name of the nuisance parameter), the second
    # is the 'effect' as a fraction, the third one is the process name. The fourth parameter
    # is optional and denotes the channl. The default '*' means that the uncertainty applies
    # to all channels in the same way.
    # Note that you can use the same name for a systematic here as for a shape
    # systematic. In this case, the same parameter will be used; shape and rate changes 
    # will be 100% correlated.

    model.add_lognormal_uncertainty('TTMSDecays_rate',      math.log(1.30), 'TTMSDecays'    )
    model.add_lognormal_uncertainty('WExclb_rate',          math.log(2.00), 'WExclb'        )
    model.add_lognormal_uncertainty('WExclc_rate',          math.log(2.00), 'WExclc'        )
    model.add_lognormal_uncertainty('WExcll_rate',          math.log(2.00), 'WExcll'        )
    model.add_lognormal_uncertainty('DY_rate',              math.log(1.30), 'DY'            )
    model.add_lognormal_uncertainty('SingleTop_rate',       math.log(1.30), 'SingleTop'     )
    model.add_lognormal_uncertainty('SingleTopW_rate',      math.log(1.30), 'SingleTopW'    )
    model.add_lognormal_uncertainty('VV_rate',              math.log(1.30), 'VV'            )
    model.add_lognormal_uncertainty('QCD_rate',             math.log(1.50), 'QCD'           )

    #model.distribution.set_distribution_parameters(benchmark, mean=0.0,width=0.0, range =[0.0,0.0])

    for p in model.processes:
        # because QCD is fully datadriven
        if p == 'QCD' : continue
	model.add_lognormal_uncertainty('lumi',        math.log(1.026), p)
                               
    return model



model = get_model()

model_summary(model)

print ("------------------------------------------------------------------")
print ("------------------------------------------------------------------")


### For max. Likelihood Fit results

print ("Run MLE")

#signal_shapes = {'': []}
signal_shapes = {benchmark:[benchmark]}
                    
fit = mle(model, input = 'data', n = 1, signal_process_groups = signal_shapes, with_covariance=False, with_error=True, ks = True, chi2 = True, options = options)

# the output is (fitted value, uncertainty)
# The first numbers in the brackets show how far we are from the nominal value (which is 0) after the fit. 
#A value of 1 would mean 1 sigma deviation. So we are below 1 sigma deviation. 
#The second numbers in the brackets illustrates the uncertainty on the fitted value, it should be below 1, 
#and a value close to 1 corresponds to "no sensitivity" on the systematic.

print ("Determine nuisance parameters and their uncertainties")

parameter_values = {}
parameter_uncert = {}

################################################
###### To quantify the changes when using ######
######   a +1 sigma nuisance parameter    ###### 
################################################
#
if doUpVariation:

    print ("Start moving one by one each nuisance parameter by 1 sigma ")

    parameter_values_systup = {}
    parameter_uncert_systup = {}

    for q in model.get_parameters([]):
        print 
        print q
        model_systup = model.copy()
	for p in model_systup.get_parameters([]):
    	    parameter_values_systup[p] = fit[benchmark][p][0][0]
    	    parameter_uncert_systup[p] = fit[benchmark][p][0][1]
            if p == q:
                # print [par_value, par_uncert] before the systUp fit
                print [q, "%.4f" %parameter_values_systup[p], "%.4f" %parameter_uncert_systup[p] ]
                if parameter_values_systup[p] < 0 : parameter_values_systup[p] = fit[benchmark][p][0][0]-fit[benchmark][p][0][1]
                else                              : parameter_values_systup[p] = fit[benchmark][p][0][0]+fit[benchmark][p][0][1]

        # fix q-parameter to its +1 sigma value
        model_systup.distribution.set_distribution_parameters(q, width = 0.0, mean = parameter_values_systup[q], range = [parameter_values_systup[q], parameter_values_systup[q]])

        # fit again but with q fixed
        fit_systup = mle(model_systup, input = 'data', n = 1, signal_process_groups = signal_shapes, with_covariance=False, with_error=True, ks = True, chi2 = True, options = options)

        # fill syst_up vectors with new values
        parameter_values_systup[q] = fit_systup[benchmark][q][0][0]
        parameter_uncert_systup[q] = fit_systup[benchmark][q][0][1]

        # print [par_value, par_uncert] after the systUp fit
        print [q, "%.4f" %parameter_values_systup[q], "%.4f" %parameter_uncert_systup[q] ]

        # to deal with the beta_signal w/o getting an error from it
        parameter_values_systup['beta_signal'] = fit[benchmark][q][0][0]

        # ouput of histos with q-systUp
        histos = evaluate_prediction(model_systup, parameter_values_systup, include_signal = True)
        write_histograms_to_rootfile(histos, 'interRegionTest/histos_postFit_SignalAndInterRegion_'+q+'.root')


# to also always get the nominal fit together with the nominal constraints
print
print "Results from a nominal fit" 
print

for p in model.get_parameters([benchmark]):
    parameter_values[p] = fit[benchmark][p][0][0]
    parameter_uncert[p] = fit[benchmark][p][0][1]

    print [p, "%.4f" %parameter_values[p], "%.4f" %parameter_uncert[p] ]

histos = evaluate_prediction(model, parameter_values, include_signal = False)
write_histograms_to_rootfile(histos, 'outputTheta_WTTnterCRSR.root')


# 2.a. Bayesian limits
# Calculate expected and observed Bayesian limits. For faster run time of this example,
# only make a few mass points. (Omitting the 'signal_procsses' parameter completely would
# process all signals defined as signal processes before; see Section "Common Parameters"
# on the theta auto intro doxygen page for details)

plot_exp, plot_obs = bayesian_limits(model, options = options, n_toy = 3000, n_data=500)

# plot_exp and plot_obs are instances of plotutil.plotdata. they contain x/y values and
# bands. You can do many things with these objects such as inspect the x/y/ban
# data, pass then to plotutil.plot routine to make pdf plots, ...
# Here, we will just create text files of the plot data. This is useful if you want
# to apply your own plotting routines or present the result in a text Table.


plot_exp.write_txt('limits_WTTinterCRSR/bayesian_limits_expected_'+benchmark+'.txt')
plot_obs.write_txt('limits_WTTinterCRSR/bayesian_limits_observed_'+benchmark+'.txt')
#plot_exp.write_txt('limits_WTTSR_Wmerged/bayesian_limits_expected_'+benchmark+'.txt')
#plot_obs.write_txt('limits_WTTSR_Wmerged/bayesian_limits_observed_'+benchmark+'.txt')

# 2.b. CLs limits
# calculate cls limit plots. The interface is very similar to bayesian_limits. However, there are a few
# more options such as the definition of the test statistic which is usually a likelihood ratio but can differ in
# which parameters are minimized and which constraints / ranges are applied during minimization.
# Here, we stay with the default which fixes beta_signal=0
# for the background only hypothesis and lets it float freely for the signal+background hypothesis.
# See cls_limits documentation for more options.

#######plot_exp, plot_obs = cls_limits(model)

# as for the bayesian limits: write the result to a text file

########plot_exp.write_txt('limits/cls_limits_expected_'+benchmark+'.txt')
########plot_obs.write_txt('limits/cls_limits_observed_'+benchmark+'.txt')

report.write_html('htmlout')

