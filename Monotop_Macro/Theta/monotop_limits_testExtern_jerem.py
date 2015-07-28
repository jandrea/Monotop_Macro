options = Options()

options.set('minimizer', 'strategy', 'newton_vanilla')

benchmark='S1Res700Inv150'

dosysttable=False

# for model building:
def get_model():
    # Read in and build the model automatically from the histograms in the root file. 
    # This model will contain all shape uncertainties given according to the templates
    # which also includes rate changes according to the alternate shapes.
    # For more info about this model and naming conventuion, see documentation
    # of build_model_from_rootfile.
    
    model = build_model_from_rootfile('thetaInOut/inputTheta_merged_AllRegions.root',include_mc_uncertainties=True)
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

    #model.distribution.set_distribution_parameters('WExclb_rate', mean=0.0,width=0.0, range =[0.0,0.0])

    for p in model.processes:
        # because QCD is fully datadriven
        if p == 'QCD' : continue
	model.add_lognormal_uncertainty('lumi',        math.log(1.026), p)
                               
    return model

model = get_model()

print ("----------------------------------------------------------------------")
print ("------------------------externalize systematics-----------------------")

#fixed_syst_list = ['match','scale', 'trig', 'fit']
fixed_syst_list = []

for fix_uncertainties in (None, fixed_syst_list ):
    print "\nuncertainties not fitted: ", fix_uncertainties
    if fix_uncertainties is None: fixed_dist = None
    else: fixed_dist = get_fixed_dist_at_values(dict([(u, 0.0) for u in fix_uncertainties]))

    print fix_uncertainties

print ("------------------------------------------------------------------")
print ("------------------------------------------------------------------")


model_summary(model)


### For max. Likelihood Fit results

print ("Run MLE")

signal_shapes = {benchmark:[benchmark]}
                    
one_sigma   = 0.6827


res    = pl_interval(model, 'data', n=1, cls = [one_sigma], signal_process_groups = signal_shapes, options = options, nuisance_constraint = fixed_dist)


print ("Determine nuisance parameters and their uncertainties")

syst_down   = 0.00000001
syst_up     = 0.00000001

if (res[benchmark][0][0]> 0.0000001):
        syst_down   = (res[benchmark][0][0]            - res[benchmark][one_sigma][0][0])/res[benchmark][0][0]
        syst_up     = (res[benchmark][one_sigma][0][1] - res[benchmark][0][0])/res[benchmark][0][0]


interval = (res[benchmark][one_sigma][0][1] - res[benchmark][one_sigma][0][0])/2
print "interval: ", interval 

#signal_fit  = signal_init_xs*res[benchmark][0][0]
#signal_down = signal_init_xs*res[benchmark][one_sigma][0][0]
#signal_up   = signal_init_xs*res[benchmark][one_sigma][0][1] 

signal_fit  = res[benchmark][0][0]
signal_down = res[benchmark][one_sigma][0][0]
signal_up   = res[benchmark][one_sigma][0][1] 


#print ["fitted cross section ", "%.5f" %signal_fit]
#print ["down variation       ", "%.5f" %signal_down]
#print ["up variation         ", "%.5f" %signal_up]

print ["up syst              ", "%.5f" %syst_up]
print ["down syst            ", "%.5f" %syst_down]

#print ["cross section ", "%.5f" %signal_fit, "+%.5f" %(signal_up-signal_fit), "-%.5f" %(signal_fit - signal_down)]


filename = "tables.tex"
texfile = open(filename,'w')

print ("-----------------------------------------------------------------------------------------")
print ("----------------------account for externalized systematics-------------------------------")

texfile.write( "the fitted beta signal   %.5f  [%.5f %.5f] \n" % ( res[benchmark][0][0], res[benchmark][one_sigma][0][0] , res[benchmark][one_sigma][0][1] ) )
texfile.write( "the fitted cross-section %.5f  [%.5f %.5f] \n" % ( signal_fit, signal_down, signal_up) )
texfile.write( "Summary of externalize uncertainties \n")
texfile.write( "\\begin{table} \n")
texfile.write( "\\begin{center} \n")
texfile.write( "\\begin{tabular}{ |c|c|c| } \n")
texfile.write( "\\hline \n")
texfile.write( "Uncertainty source & up (\\%) & down (\\%) \\\\ \n")
texfile.write( "\\hline \n")

total_error_up   = (signal_up - signal_fit)**2
total_error_down = (signal_fit - signal_down)**2

fixed_uncert_error_up = 0
fixed_uncert_error_down = 0


for p in fix_uncertainties:

        #print [p, math.sqrt(total_error_up), math.sqrt(total_error_down)]
        model_tmp = model.copy()
        p_mean = model_tmp.distribution.get_distribution(p)['mean']
        p_width = model_tmp.distribution.get_distribution(p)['width']
        #-------------------------
        #for up shift
        #-------------------------
        nuisance_prior_toys = get_fixed_dist_at_values({p: p_mean + p_width})   
        res_up = mle(model_tmp, 'toys-asimov:0.77363', 1, nuisance_prior_toys = nuisance_prior_toys, nuisance_constraint = fixed_dist, options=options)
        beta_signal_fitted_up = res_up[benchmark]['beta_signal'][0][0]
        shift_plus = beta_signal_fitted_up - res[benchmark][0][0]
        #xs_shift_plus= signal_init_xs*shift_plus

        #-------------------------
        #for down shift
        #-------------------------
        nuisance_prior_toys = get_fixed_dist_at_values({p: p_mean - p_width})
        res_down = mle(model, 'toys-asimov:0.77363', 1, nuisance_prior_toys = nuisance_prior_toys, nuisance_constraint = fixed_dist, options=options)
        beta_signal_fitted_down = res_down[benchmark]['beta_signal'][0][0]
        shift_minus    = beta_signal_fitted_down - res[benchmark][0][0] 
        #xs_shift_minus = signal_init_xs*shift_minus
                                        
        print [p, beta_signal_fitted_up, beta_signal_fitted_down]
        #print [p, shift_plus, shift_minus]
            
                
        #texfile.write( "%s & %.5f  & %.5f \\\\ \n" % (p, shift_plus*100, shift_minus*100) )
        #total_error_up   += xs_shift_plus**2
        #total_error_down += xs_shift_minus**2
        #fixed_uncert_error_up   += xs_shift_plus**2
        #fixed_uncert_error_down += xs_shift_minus**2
                
        texfile.write( "%s & %.5f  & %.5f \\\\ \n" % (p, shift_plus*100, shift_minus*100) )
        total_error_up   += shift_plus**2
        total_error_down += shift_minus**2
        fixed_uncert_error_up   += shift_plus**2
        fixed_uncert_error_down += shift_minus**2

print [ math.sqrt(total_error_up), math.sqrt(total_error_down)] 
texfile.write( "\\end{tabular}\n")
texfile.write( "\\end{center} \n")
texfile.write( "\\end{table} \n")

total_error_up   =  math.sqrt(total_error_up)
total_error_down =  math.sqrt(total_error_down)

print ["fitted cross section with externalized systematics", "%.5f" %signal_fit, "+ %.5f" %total_error_up, "- %.5f" %total_error_down]


texfile.write( "fitted cross section with externalized systematics %.4f+%.4f-%.4f pb \n" % (signal_fit, total_error_up, total_error_down) )
texfile.write( "up variation   %.5f \n" % total_error_up)


print ("------------------------------------------------------------------")
print ("------------------------------------------------------------------")
print "Signal fit: ", signal_fit 


### For max. Likelihood Fit results

print ["start fit mle"]
fit = mle(model, input = 'data', n = 1, signal_process_groups = signal_shapes, with_covariance=False, with_error=True, ks = True, chi2 = True, options = options)
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

parameter_values['beta_signal'] =  res[benchmark][0][0]

histos = evaluate_prediction(model, parameter_values, include_signal = False)
#write_histograms_to_rootfile(histos, 'thetaInOut/outputTheta_merged_CRsOnly.root')
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




if dosysttable:


    texfile.write( " \n")
    texfile.write( " \n")
    texfile.write( " \n")

    texfile.write( "\\begin{table} \n")
    texfile.write( "\\begin{center} \n")
    texfile.write( "\\begin{tabular}{ |c|c|c| } \n")
    texfile.write( "\\hline \n")
    texfile.write( "Uncertainty source & pb & \\%  \\\\ \n")
    texfile.write( "\\hline \n")

    print ("--------------------------------------------------------------------")
    print ("-----------------------do systematic table--------------------------")
    syst = 0
    tot_uncert =0

    print ["--------------------------------"]
    print ("Determine the impact of each systematic")
    for p in model.get_parameters([]):
        if(p == 'QCD_rate'): continue
        if(p == 'WExcll_rate'): continue
        if(p == 'WExclb_rate'): continue
        if(p == 'WExclc_rate'): continue
        if(p == 'TTMSDecays_rate'): continue
        if(p == 'DY_rate'): continue         
        if(p == 'SingleTop_rate'): continue  
        if(p == 'SingleTopW_rate'): continue 
        if(p == 'VV_rate'): continue         
        #if(p == 'trig'): continue
        if(p == 'toppt'): continue

        excluded_syst=False
        for q in fixed_syst_list:
            if(p==q): excluded_syst=True
        if(excluded_syst==True): continue
        print p
        model_syst = model.copy()
        model_syst.distribution.set_distribution_parameters(p, width = 0.0, mean = parameter_values[p], range = [parameter_values[p], parameter_values[p]])
        res_syst = pl_interval(model_syst, 'data', n=1, cls = [one_sigma], signal_process_groups = signal_shapes, options = options  )


        interval_syst = (res_syst[benchmark][one_sigma][0][1] - res_syst[benchmark][one_sigma][0][0])/2
        syst  = (abs(interval**2 - interval_syst**2))**(0.5)
        tot_uncert = tot_uncert + (syst)**2
        #print ["syst effect of (%)  ", p, "%.2f" %(syst*100)]
        #print ["syst effect of (pb) ", p, "%.2f" %(syst*signal_fit)]
        print [p, "%.4f" %(syst*signal_fit), "%.4f" %(syst*100)]

        texfile.write( "%s & %.5f  & %.5f \\\\ \n" % (p, syst*signal_fit, syst*100) )

            #print ["--------------------------------"]

    #print ["total syst down/up" ,"%.4f" % total_down**(0.5), "%.4f" %total_up**(0.5)]
    tot_uncert = tot_uncert**(0.5)
    print ["total uncert %" , tot_uncert*100]
    #tot_uncert = tot_uncert*signal_fit
    print ["total uncert pb" , tot_uncert*signal_fit]
    texfile.write( "\\hline \n")
    texfile.write( "total syst & %.5f  & %.5f \\\\ \n" % (tot_uncert*signal_fit, tot_uncert*100 ) )
    texfile.write( "\\hline \n")

    texfile.write( "\\end{tabular}\n")
    texfile.write( "\\end{center} \n")
    texfile.write( "\\end{table} \n")
    texfile.write( " \n")
    texfile.write( " \n")
    syst_error = ( (tot_uncert*signal_fit)**(2) + fixed_uncert_error_up**(2))**(0.5)
    texfile.write( " %.5f \n" % total_error_up)
    texfile.write( " %.5f \n" % syst_error)
    texfile.write( " %.5f \n" % (total_error_up**(2) - syst_error**(2))**(0.5) )
    stat_error = 0
    texfile.write( "final cross section is %.2f \n" % signal_fit )
    texfile.write( "total uncertainty %.1f %.2f\n" % ( total_error_up, total_error_down ) )
    texfile.write( "stat %.2f \n" % (total_error_up**(2) - syst_error**(2))**(0.5) )
    texfile.write( "syst %.2f \n" % (syst_error) ) 


#print ("------------------------------------------------------------------")
#print ("------------------------------------------------------------------")


## pour propager dans l'externalisation calcul de limite
#plot_exp_bayes, plot_obs_bayes = bayesian_limits(model, 'all', n_toy = 2000, n_data = 20,nuisance_constraint=fixed_dist)
#plot_exp_bayes.write_txt('test_fix_syst/bayesian_limits_expected.txt')
#plot_obs_bayes.write_txt('test_fix_syst/bayesian_limits_observed.txt')



# 2.a. Bayesian limits
# Calculate expected and observed Bayesian limits. For faster run time of this example,
# only make a few mass points. (Omitting the 'signal_procsses' parameter completely would
# process all signals defined as signal processes before; see Section "Common Parameters"
# on the theta auto intro doxygen page for details)

###plot_exp, plot_obs = bayesian_limits(model, options = options, n_toy = 2000, n_data=200)

# plot_exp and plot_obs are instances of plotutil.plotdata. they contain x/y values and
# bands. You can do many things with these objects such as inspect the x/y/ban
# data, pass then to plotutil.plot routine to make pdf plots, ...
# Here, we will just create text files of the plot data. This is useful if you want
# to apply your own plotting routines or present the result in a text Table.


#plot_exp.write_txt('limits_CutnCount_ATLASSel/bayesian_final_limits_expected_'+benchmark+'.txt')
#plot_obs.write_txt('limits_CutnCount_ATLASSel/bayesian_final_limits_observed_'+benchmark+'.txt')
###plot_exp.write_txt('limits_ShapeAnalysis_AllRegions/bayesian_final_limits_expected_'+benchmark+'.txt')
###plot_obs.write_txt('limits_ShapeAnalysis_AllRegions/bayesian_final_limits_observed_'+benchmark+'.txt')

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

