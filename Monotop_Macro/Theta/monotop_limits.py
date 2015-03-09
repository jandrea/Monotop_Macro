options = Options()

options.set('minimizer', 'strategy', 'newton_vanilla')

benchmark='S1res1300inv100'
#benchmark='S1res1000inv100'
#benchmark='S1res1000inv800'
#benchmark='S4inv400'
#benchmark='S4inv600'
#benchmark='S4inv700'
#benchmarks='S*inv*'

# for model building:
def get_model():
    # Read in and build the model automatically from the histograms in the root file. 
    # This model will contain all shape uncertainties given according to the templates
    # which also includes rate changes according to the alternate shapes.
    # For more info about this model and naming conventuion, see documentation
    # of build_model_from_rootfile.
    
    model = build_model_from_rootfile('inputTheta_merged.root',include_mc_uncertainties=True)
#    model = build_model_from_rootfile('inputTheta_signalregion.root',include_mc_uncertainties=True)
#    model = build_model_from_rootfile('inputTheta_Wregion.root',include_mc_uncertainties=True)
    
    # If the prediction histogram is zero, but data is non-zero, teh negative log-likelihood
    # is infinity which causes problems for some methods. Therefore, we set all histogram
    # bin entries to a small, but positive value:
    model.fill_histogram_zerobins()

    # define what the signal processes are. All other processes are assumed to make up the 
    # 'background-only' model.
    model.set_signal_processes('S1*')



    # Add some lognormal rate uncertainties. The first parameter is the name of the
    # uncertainty (which will also be the name of the nuisance parameter), the second
    # is the 'effect' as a fraction, the third one is the process name. The fourth parameter
    # is optional and denotes the channl. The default '*' means that the uncertainty applies
    # to all channels in the same way.
    # Note that you can use the same name for a systematic here as for a shape
    # systematic. In this case, the same parameter will be used; shape and rate changes 
    # will be 100% correlated.

    model.add_lognormal_uncertainty('TTbarMadgraph_rate',   math.log(1.30), 'TTbarMadgraph' )
    model.add_lognormal_uncertainty('WExclb_rate',          math.log(2.00), 'WExclb'        )
    model.add_lognormal_uncertainty('WExclc_rate',          math.log(2.00), 'WExclc'        )
    model.add_lognormal_uncertainty('WExcll_rate',          math.log(2.00), 'WExcll'        )
    model.add_lognormal_uncertainty('DY10To50_rate',        math.log(1.30), 'DY10To50'      )
    model.add_lognormal_uncertainty('DY50_rate',            math.log(1.30), 'DY50'          )
    model.add_lognormal_uncertainty('Ts_rate',              math.log(1.30), 'Ts'            )
    model.add_lognormal_uncertainty('Tt_rate',              math.log(1.30), 'Tt'            )
    model.add_lognormal_uncertainty('TtW_rate',             math.log(1.30), 'TtW'           )
    model.add_lognormal_uncertainty('Tbart_rate',           math.log(1.30), 'Tbart'         )
    model.add_lognormal_uncertainty('TbartW_rate',          math.log(1.30), 'TbartW'        )
    model.add_lognormal_uncertainty('WZ_rate',              math.log(1.30), 'WZ'            )
    model.add_lognormal_uncertainty('ZZ_rate',              math.log(1.30), 'ZZ'            )
    model.add_lognormal_uncertainty('WW_rate',              math.log(1.30), 'WW'            )

    #model.distribution.set_distribution_parameters('WExclb_rate', mean=0.0,width=0.0, range =[0.0,0.0])

    for p in model.processes:
	if p == 'QCD' : continue
	if p == 'QCDA': continue
	if p == 'QCDB': continue
	if p == 'QCDC': continue
	if p == 'QCDD': continue
	model.add_lognormal_uncertainty('lumi',        math.log(1.026), p)

    for q in model.processes:
        if q == 'QCD' : model.add_lognormal_uncertainty('QCDA_rate',            math.log(1.50), 'QCD'           )
        if q == 'QCDA': model.add_lognormal_uncertainty('QCDA_rate',            math.log(1.50), 'QCDA'          )
        if q == 'QCDB': model.add_lognormal_uncertainty('QCDB_rate',            math.log(1.50), 'QCDB'          )
        if q == 'QCDC': model.add_lognormal_uncertainty('QCDC_rate',            math.log(1.50), 'QCDC'          )
        if q == 'QCDD': model.add_lognormal_uncertainty('QCDD_rate',            math.log(1.50), 'QCDD'          )
                               
    return model



model = get_model()

model_summary(model)

print ("------------------------------------------------------------------")
print ("------------------------------------------------------------------")


### For max. Likelihood Fit results

print ("Run MLE")

#signal_shapes = {'S4inv400':['S4inv400'],'S4inv600':['S4inv600'],'S4inv700':['S4inv700']}
signal_shapes = {benchmark:[benchmark]}
#signal_shapes = {'S1res500inv100':['S1res500inv100']}
#signal_shapes = {'S4inv400':['S4inv400']}
#signal_shapes = model.signal_process_groups
                    
fit = mle(model, input = 'data', n = 1, signal_process_groups = signal_shapes, with_covariance=False, with_error=True, ks = True, chi2 = True, options = options)

# the output is (fitted value, uncertainty)
# The first numbers in the brackets show how far we are from the nominal value (which is 0) after the fit. 
#A value of 1 would mean 1 sigma deviation. So we are below 1 sigma deviation. 
#The second numbers in the brackets illustrates the uncertainty on the fitted value, it should be below 1, 
#and a value close to 1 corresponds to "no sensitivity" on the systematic.

print ("Determine nuisance parameters and their uncertainties")

parameter_values = {}
parameter_uncert = {}

for p in model.get_parameters([benchmark]):
    parameter_values[p] = fit[benchmark][p][0][0]
    parameter_uncert[p] = fit[benchmark][p][0][1]
    print [p, "%.4f" %parameter_values[p], "%.4f" %parameter_uncert[p] ]

print ("Create postfit histograms")

histos = evaluate_prediction(model, parameter_values, include_signal = False)
write_histograms_to_rootfile(histos, 'histos-mle.root')


# 2.a. Bayesian limits
# Calculate expected and observed Bayesian limits. For faster run time of this example,
# only make a few mass points. (Omitting the 'signal_procsses' parameter completely would
# process all signals defined as signal processes before; see Section "Common Parameters"
# on the theta auto intro doxygen page for details)

plot_exp, plot_obs = bayesian_limits(model, options = options, n_toy = 1000, n_data=200)

# plot_exp and plot_obs are instances of plotutil.plotdata. they contain x/y values and
# bands. You can do many things with these objects such as inspect the x/y/ban
# data, pass then to plotutil.plot routine to make pdf plots, ...
# Here, we will just create text files of the plot data. This is useful if you want
# to apply your own plotting routines or present the result in a text Table.


plot_exp.write_txt('limits/bayesian_limits_expected_.txt')
plot_obs.write_txt('limits/bayesian_limits_observed_.txt')

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

