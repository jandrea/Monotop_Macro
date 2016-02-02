options = Options()

options.set('minimizer', 'strategy', 'newton_vanilla')

benchmark='S1Res900Inv100'

doextern=False
dobreakdown=True

# for model building:
def get_model():
    # Read in and build the model automatically from the histograms in the root file. 
    # This model will contain all shape uncertainties given according to the templates
    # which also includes rate changes according to the alternate shapes.
    # For more info about this model and naming conventuion, see documentation
    # of build_model_from_rootfile.
    
    model = build_model_from_rootfile('inputTheta_merged_AllRegions_RES.root',include_mc_uncertainties=True)
    #model = build_model_from_rootfile('thetaInOut/inputTheta_merged_CRsOnly.root',include_mc_uncertainties=True)
    
    # If the prediction histogram is zero, but data is non-zero, teh negative log-likelihood
    # is infinity which causes problems for some methods. Therefore, we set all histogram
    # bin entries to a small, but positive value:
    model.fill_histogram_zerobins(0.0000000001)

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

    model.add_lognormal_uncertainty('TTMSDecays_rate',      math.log(1.30), 'TTMSDecays'    )
    model.add_lognormal_uncertainty('WExclb_rate',          math.log(2.00), 'WExclb'        )
    model.add_lognormal_uncertainty('WExclc_rate',          math.log(2.00), 'WExclc'        )
    model.add_lognormal_uncertainty('WExcll_rate',          math.log(2.00), 'WExcll'        )
    model.add_lognormal_uncertainty('DY_rate',              math.log(1.30), 'DY'            )
    model.add_lognormal_uncertainty('SingleTop_rate',       math.log(1.30), 'SingleTop'     )
    model.add_lognormal_uncertainty('SingleTopW_rate',      math.log(1.30), 'SingleTopW'    )
    model.add_lognormal_uncertainty('VV_rate',              math.log(1.30), 'VV'            )
    model.add_lognormal_uncertainty('QCD_rate',             math.log(1.50), 'QCD'           )
    #model.add_lognormal_uncertainty('S4Inv600_rate',       math.log(1.10), 'S4Inv600'      )

    #model.distribution.set_distribution_parameters('WExclb_rate', mean=0.0,width=0.0, range =[0.0,0.0])

    for p in model.processes:
        # because QCD is fully datadriven
        if p == 'QCD' : continue
	model.add_lognormal_uncertainty('lumi',        math.log(1.026), p)
                               
    return model

model = get_model()



if doextern:

    print ("----------------------------------------------------------------------")
    print ("------------------------externalize systematics-----------------------")

    #fixed_syst_list = ['match','scale', 'trig', 'fit']
    fixed_syst_list = ['matching', 'scale', 'PDF', 'toppt']

    for fix_uncertainties in (None, fixed_syst_list ):
        print "\nuncertainties not fitted: ", fix_uncertainties
        if fix_uncertainties is None: fixed_dist = None
        else: fixed_dist = get_fixed_dist_at_values(dict([(u, 0.0) for u in fix_uncertainties]))

        print fix_uncertainties

    print ("----------------------------------------------------------------------")
    print ("----------------------------------------------------------------------")

model_summary(model)


### For max. Likelihood Fit results

print ("Run MLE")

signal_shapes = {benchmark:[benchmark]}

### For max. Likelihood Fit results

print ["start fit mle"]
if doextern: fit = mle(model, input = 'data', n = 1, signal_process_groups = signal_shapes, with_covariance=False, with_error=True, ks = True, chi2 = True, options = options, nuisance_constraint = fixed_dist)
else: fit = mle(model, input = 'data', n = 1, signal_process_groups = signal_shapes, with_covariance=False, with_error=True, ks = True, chi2 = True, options = options)


# the output is (fitted value, uncertainty)
# The first numbers in the brackets show how far we are from the nominal value (which is 0) after the fit. 
#A value of 1 would mean 1 sigma deviation. So we are below 1 sigma deviation. 
#The second numbers in the brackets illustrates the uncertainty on the fitted value, it should be below 1, 
#and a value close to 1 corresponds to "no sensitivity" on the systematic.
print ["mle fit done"]

print ("Determine nuisance parameters and their uncertainties")
parameter_values = {}
parameter_uncert = {}


for p in model.get_parameters([benchmark]):
    parameter_values[p] = fit[benchmark][p][0][0]
    parameter_uncert[p] = fit[benchmark][p][0][1]

    print [p, "%.4f" %parameter_values[p], "%.4f" %parameter_uncert[p] ]

histos = evaluate_prediction(model, parameter_values, include_signal = False)
write_histograms_to_rootfile(histos, 'outputTheta_AllRegions.root')
#write_histograms_to_rootfile(histos, 'outputTheta_Extern_AllRegions.root')
#write_histograms_to_rootfile(histos, 'thetaInOut/outputTheta_merged_'+benchmark+'_AllRegions.root')
#write_histograms_to_rootfile(histos, 'thetaInOut/outputTheta_testExtern_merged_'+benchmark+'_AllRegions.root')

print ("------------------------------------------------------------------")
print ("------------------------------------------------------------------")


################################################
###### To quantify the changes when using ######
######   a +1 sigma nuisance parameter    ###### 
################################################
#
#for q in model.get_parameters([benchmark]):
#    for p in model.get_parameters([benchmark]):
#        parameter_values[p] = fit[benchmark][p][0][0]
#        parameter_uncert[p] = fit[benchmark][p][0][1]
#
#        if p == q : 
#            parameter_values[p] = fit[benchmark][p][0][0]+fit[benchmark][p][0][1]
#            print [p, "%.4f" %parameter_values[p], "%.4f" %parameter_uncert[p] ]
#
#    histos = evaluate_prediction(model, parameter_values, include_signal = False)
#    #write_histograms_to_rootfile(histos, 'thetaInOut/histos_postFit_AllRegions_'+q+'.root')
#    write_histograms_to_rootfile(histos, 'thetaInOut/histos_postFit_CRsOnly_'+q+'.root')
#
#
#################################################

if doextern: plot_exp_bayes, plot_obs_bayes = bayesian_limits(model, options = options, n_toy = 10000, n_data = 2500, nuisance_constraint = fixed_dist)
else: plot_exp_bayes, plot_obs_bayes = bayesian_limits(model, options = options, n_toy = 3000, n_data = 500)
#plot_exp_bayes.write_txt('bayesian_limits_expected_nominal.txt')
#plot_obs_bayes.write_txt('bayesian_limits_observed_nominal.txt')
plot_exp_bayes.write_txt('limitsBreakdown/bayesian_limits_expected_'+benchmark+'.txt')
plot_obs_bayes.write_txt('limitsBreakdown/bayesian_limits_observed_'+benchmark+'.txt')

if dobreakdown:

    print ("------------------------------------------------------------------------")
    print ("--------------------- do breakdown of systematics ----------------------")
    print ("--------------- determine the impact of each systematics ---------------")
    print ("------------------------------------------------------------------------")

    for p in model.get_parameters([benchmark]):
        #if(p != 'S1Res700Inv150_rate'): continue
        #if(p != 'PDF'): continue
        #if(p == 'QCD_rate'): continue
        #if(p == 'WExcll_rate'): continue
        #if(p == 'WExclb_rate'): continue
        #if(p == 'WExclc_rate'): continue
        #if(p == 'TTMSDecays_rate'): continue
        #if(p == 'DY_rate'): continue         
        #if(p == 'SingleTop_rate'): continue  
        #if(p == 'SingleTopW_rate'): continue 
        #if(p == 'VV_rate'): continue         
        if(p == 'beta_signal'): continue

        excluded_syst=False
        for q in model.get_parameters([benchmark]):
            if(p==q): excluded_syst=True
        if(excluded_syst==False): continue # i.e. if it is a background rate and not a systematics
        print p
        model_syst = model.copy()
        model_syst.distribution.set_distribution_parameters(p, width = 0.0, mean = parameter_values[p], range = [parameter_values[p], parameter_values[p]])


        ## pour propager dans l'externalisation calcul de limite
        plot_exp_bayes_nosyst, plot_obs_bayes_nosyst = bayesian_limits(model_syst, options = options, n_toy = 3000, n_data = 500)
        plot_exp_bayes_nosyst.write_txt('limitsBreakdown/bayesian_limits_expected_'+benchmark+'_wo_'+p+'.txt')
        plot_obs_bayes_nosyst.write_txt('limitsBreakdown/bayesian_limits_observed_'+benchmark+'_wo_'+p+'.txt')


report.write_html('htmlout')

