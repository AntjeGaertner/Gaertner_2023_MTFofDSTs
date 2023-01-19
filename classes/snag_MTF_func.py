#!/usr/bin/python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.optimize as opt
from scipy.interpolate import interp1d
import statsmodels.formula.api as smf
from statsmodels.tools.eval_measures import rmse as rmse_smf
from statsmodels.stats.stattools import durbin_watson as durbin_watson_test
from statsmodels.stats.diagnostic import het_white
import sys
import os.path

#sys.path.append('/Users/antje/Update_PIK/')

from classes import boundaries as bd
from classes import functions as fcy

plot_setup = bd.plot_setup()
params = {'mathtext.default':'regular'}
plt.rcParams.update(params)

# Coefficients (multiplication!!!)
cm_to_m    = 0.01
m_to_cm    = 100
m2_to_ha   = 0.0001
acre_to_ha = 0.404686
inch_to_cm = 2.54
feet_to_m  = 0.3048


def MTF_snags(x, y, max_x, k=np.nan, plot=False, print_out=False):
     
    """
    x needs to reflect the middle between data recording intervals.
    """
   
    np.seterr(divide = 'ignore')
    # Compute relative share of snags remaining/Users/antje/Downloads/wpd_datasets (1).csv
    snags_remain_endt = 1-np.sum(y[:-1]-y[1:])
    if print_out:
        print('Snags remaining: ', "%.2f"%snags_remain_endt, '%')

    if snags_remain_endt > 0.01:
        if print_out:
            print('Computing t_99.')
        t_99 = -np.log(0.01)/k
        if print_out:
            print('t_99: ',t_99)
        
        #-- interpolate the data to reduce noise
#        x_intp = np.linspace(0,x[-1], 1000)
#        f_y = interp1d(x=x, y=y, kind='linear', fill_value='extrapolate')
        
        # replace the original input
#        y = f_y(x_intp)
#        x = x_intp
        
        #-- T1 from k
        T1 = 1/k
        if print_out:
            print('T1:', "%.2f"%T1)
        
        #-- T2
        yy = y[:-1]-y[1:] # Differences based on data only
        
        if plot:
            plt.bar(x[1:], yy, color='C0')
            plt.plot(x[1:], yy, 'r--')
            plt.title('Differences between timesteps')
            plt.ylabel('t-t+1')
            plt.xlabel('time (years)')
            plt.show()

        T2=np.nansum(yy*x[1:])/np.nansum(yy)
        if print_out:
            print('T2:', "%.2f"%T2)

        #-- T3
        timestep = 1000
        xq = np.linspace(0, x[-1], timestep)

        diff_interp = np.interp(x=xq, xp=x[1:], fp=yy)

        T3 = np.nansum(diff_interp*xq)/np.nansum(diff_interp)
        if print_out:
            print('T3:', "%.2f"%T3)
        
        #-- T4
        # Now we need t_99
        # Computing new x values based on this
        t_99 = 100
        xq = np.linspace(0, t_99, timestep)
        f = interp1d(x=x[1:], y=yy, kind='linear', fill_value='extrapolate')
        diff_interp = f(xq)

        T4 = np.nansum(diff_interp[diff_interp>0]*xq[diff_interp>0])/np.nansum(diff_interp[diff_interp>0])
        if print_out:
            print('T4:', "%.2f"%T4)
        
        
        #-- T5 - log transforming the values
        # Dropping negative change values to perform log, cause there shouldn't be suddenly more snags
        f_log = interp1d(x=x[1:][yy>=0], y=np.log(yy[yy>=0]), kind='linear', fill_value='extrapolate')
        
        diff_interp = np.exp(f_log(xq))

        T5 = np.nansum(diff_interp*xq)/np.nansum(diff_interp)
        if print_out:
            print('T5:', "%.2f"%T5)
        if plot:
            plt.bar(x[1:], yy, color='C0')
            plt.plot(xq, diff_interp, 'r--')
            plt.show()

    elif snags_remain_endt < 0.01:
	# This is relevant for Vanderwel 2006 or Storaunet 2002
        if print_out:
            print('No need to compute t_99.')

        #-- T1 from k
        T1 = 1/k
        if print_out:
            print('T1:', "%.2f"%T1)

        #-- T2 
        yy = y[:-1]-y[1:] # Differences based on data only
        if plot:
            plt.bar(x[1:], yy, color='C0')
            plt.plot(x[1:], yy, 'r--')
            plt.title('Differences between timesteps')
            plt.ylabel('t-t+1')
            plt.xlabel('time (years)')
            plt.show()

        T2=np.nansum(yy*x[1:])/np.nansum(yy)
        if print_out:
            print('T2:', "%.2f"%T2)

        #-- T3
        timestep = 100
        xq = np.linspace(0, x[-1], timestep)        

        diff_interp = np.interp(x=xq, xp=x[1:], fp=yy)

        T3 = sum(diff_interp*xq)/sum(diff_interp)
        if print_out:
            print('T3:', "%.2f"%T3)


        #-- T4
        xq = np.linspace(0, x[-1], timestep)        

        diff_interp = np.interp(x=xq, xp=x[1:], fp=yy)

        T4 = sum(diff_interp*xq)/sum(diff_interp)
        if print_out:
            print('T4:', "%.2f"%T4)

        #-- T5
        diff_interp = np.exp(np.interp(x=xq, xp=x[1:], fp=np.log(yy)))

        T5 = np.nansum(diff_interp*xq)/np.nansum(diff_interp)
        if print_out:
            print('T5:', "%.2f"%T5)
        if plot:
            plt.bar(x[1:], yy, color='C0')
            plt.plot(xq, diff_interp, 'r--')
            plt.show()
        np.seterr(divide = 'warn')
    return [T1, T2, T3, T4, T5]

    
def fit_neagtive_SEM(x_vals, y_vals):

    def func(x, k):
        return (np.exp(-k*x))
    
    #-- Fit single negativ exponential model to the data
    x = x_vals # numpy array!
    y = y_vals # numpy array!

    popt, pcov = opt.curve_fit(func, x, y,
                               maxfev= 100000
                              )
    residuals = y - func(x, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y -np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)

    print('R2=',np.round(r_squared,2), *popt)

    Y = func(x, *popt)

    k = popt[0] # Store k
    
    return [Y, k]
    
def fit_function(func, x_vals, y_vals):

    #-- Fit model to the data
    x = x_vals # numpy array!
    y = y_vals # numpy array!

    popt, pcov = opt.curve_fit(func, x, y,
                               maxfev= 100000
                              )
    residuals = y - func(x, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y -np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)

    print('R2=',np.round(r_squared,2), *popt)

    Y = func(x, *popt)

    return [Y, popt, r_squared]
    

def weight_by_single(vals, weight):
    try:
        shape = np.shape(vals)[1]
        weight_ts = np.repeat(weight, shape).reshape((len(vals),shape))
        weighted_ts = np.sum(vals * weight_ts,axis=0)/np.sum(weight)
    except:
        weight_ts = weight
        weighted_ts = np.sum(vals * weight_ts)/np.sum(weight)
    return weighted_ts


def snag_input_wrapper(data, iteratable, interp = True, fig_name = np.nan, MTF_method = 1, 
                       reference=np.nan):
    """
    Takes a dataframe and performs linear interpolation of the data.
    
    df         --> Pandas Datafreme with multiindex columns of level 0 = <iteratable(s)>, level 1 = [X, Y]
    iteratable --> list of species names or DBH classes that define <iteratable(s)> in multiindex column 
                   level 0of df, ndim = 1;
    fig_name   --> additiona name to figure title such as the figure number of original article, defaults to NaN
    MTF_method --> Method accepted as result from MTF_snags; slice(None,None) - Shows the result matrixes in full
    reference  --> Reference of the article the data is taken from
    """
    # Find maximum value of the data and contruct endpoint for the interpolation
    val = np.nanmax(data.values)
    bins = np.arange(0,200,50)
    bins[0] = 25
    end = bins[bins>val][0]
    
    # Values for the interpolation
    xs = np.arange(0,end,0.1)
    
    # Pre-allocate result matrices / lists

    result_tree = np.empty((len(iteratable),5))+np.nan
    Y = np.empty((len(iteratable),len(xs)))
    Y_original_len = []
    X_original_len = []
    f_Y = []
    
    # Loop over the interratable(s)
    idx=0
    jdx=0
    plt.figure(figsize=(15,8))
    for i in range(len(iteratable)):
        sp = iteratable[i]

        x = data.iloc[:,idx].dropna().values
        y = data.iloc[:,idx+1].dropna().values

        # Store original values
        Y_original_len.append(y)
        X_original_len.append(x)

        f_Y.append(interp1d(x,y, kind='linear',fill_value='extrapolate'))
        Y[jdx,:] = f_Y[jdx](xs)
        Y[jdx][Y[jdx]<0] = 0
        plt.plot(xs, Y[jdx], label=sp)
        plt.scatter(x,y, label=sp)

        idx+=2
        jdx+=1

    plt.title(reference+': '+fig_name)
    plt.legend()
    plt.xlabel('Time since death [years]')
    plt.ylabel('Proportion of snags \n remaining')
    plt.show()    

    jdx=0
    for i in range(len(iteratable)):
        print(iteratable[i])
        if interp == True:
            result_tree[i,MTF_method] = MTF_snags(xs, Y[jdx],np.nan)[MTF_method]
        else:
            result_tree[i,MTF_method] = MTF_snags(X_original_len[i][1:], Y_original_len[i][:-1],np.nan)[MTF_method]
        jdx+=1
    
    return xs, Y, X_original_len, Y_original_len, result_tree



def save_MTF_input(x, y,
                   species,
                   species_scientific,
                   dt, reference,
                   y_unit = 'count',
                   species_ntree = np.nan,
                   DBH_class = np.nan,
                   DBH_min = np.nan,
                   DBH_max = np.nan,
                   DBH_mean = np.nan,
                   DBH_ntree = np.nan,
                   cmass_whole = np.nan,
                   cmass_ts = np.nan,
                   site_name = np.nan):
    """
    Saves input data from MTF

    if species > 1 | DBH_class > 1, then x, y, species, DBH_class need to be lists

    Site information in site.csv to be found with reference
    """
    # x,y,species, DBH, need to be lists if n>1 for species & DBH class

    cols = ['x', 'y', 'y.unit', 'species', 'species_scientific', 'species_ntree', 'DBH_class', 'DBH_min', 'DBH_max',
            'DBH_mean', 'DBH_ntree', 'cmass_whole', 'cmass_ts', 'DT_method', 'site_name', 'Reference']

    print(all(type(i) == list for i in [x, y, species, DBH_mean]))
    print('species',type(species))
    print('DBH_mean',type(DBH_mean))
    print('x',type(x))
    print('y',type(y))

    if all(type(i) == list for i in [x, y, species]) & (type(y_unit) != list):
        df_list = []
        #print('site_name[i]',type(site_name[0]))
        for i in range(len(species)):
            df_list.append(pd.DataFrame({'x':x[i], 'y':y[i], 'y.unit':y_unit,
                                         's':species[i], 'ss':species_scientific[i],
                                         'species_ntree':species_ntree[i],
                                         'DBH_class': DBH_class[i], 'DBH_min':DBH_min[i], 'DBH_max':DBH_max[i],
                                         'DBH_mean':DBH_mean[i], 'DBH_ntree':DBH_ntree[i],
                                         'cmass_whole':cmass_whole[i], 'cmass_ts':cmass_ts[i],
                                         'dt':dt, 'site_name':site_name[i], 'Reference':reference
                                         
                                        }))
            df = pd.concat(df_list)
            df.columns = cols
            
    elif type(y_unit) == list:
        df_list = []
        for i in range(len(y_unit)):
                df_list.append(pd.DataFrame({'x':x[i], 'y':y[i], 'y.unit':y_unit[i],
                                             's':species[i], 'ss':species_scientific[i],
                                             'species_ntree':species_ntree[i],
                                             'DBH_class': DBH_class[i], 'DBH_min':DBH_min[i], 'DBH_max':DBH_max[i],
                                             'DBH_mean':DBH_mean[i], 'DBH_ntree':DBH_ntree[i],
                                             'cmass_whole':cmass_whole[i], 'cmass_ts':cmass_ts[i],
                                             'dt':dt, 'site_name':site_name[i], 'Reference':reference
                                             
                                            }))
                df = pd.concat(df_list)
                df.columns = cols
                
    else:
        df = pd.DataFrame({'x':x, 'y':y, 'y.unit':y_unit,
                           's':species, 'ss':species_scientific, 'species_ntree':species_ntree, 'DBH_class':DBH_class,
                           'DBH_min':DBH_min, 'DBH_max':DBH_max, 'DBH_mean':DBH_mean,
                            'DBH_ntree':DBH_ntree, 'cmass_whole':cmass_whole, 'cmass_ts':cmass_ts,
                           'dt':dt, 'site_name':site_name, 'Reference':reference})
        
        df.columns = cols
        
    if os.path.isfile('210910_MTF_input.csv') == True:
        file = pd.read_csv('210910_MTF_input.csv', index_col=None, encoding='utf-8-sig')
        if file.Reference.str.contains(reference).any():
            file = file[~file.Reference.str.contains(reference)]
            
        file.reset_index(drop=True, inplace=True)
        df.reset_index(drop=True, inplace=True)
        df_new = pd.concat([file, df], axis=0, sort=False)
        df_new = df_new[cols]
        
        # Store
        df_new.to_csv('210910_MTF_input.csv', index=False, encoding='utf-8-sig')
    return df_new
    
    

def MTFcount2MTFcmass(MTF_count): # SITE LEVEL!!!!
    # These coefficents are from the regression between ntree MTF and ntree BA where both where available n=10
    
    intercept  = 0
    slope_site_count2cmass = 1.9462            
    
    MTF_cmass = intercept+MTF_count*slope_site_count2cmass
    
    return MTF_cmass
    
def MTFcount2MTFcmass_species(MTF_count): # SPECIES LEVEL!!!
    # These coefficents are from the regression between ntree MTF and ntree BA where both where available n=10

    intercept = 0
    slope_species_count2cmass = 1.4802            

    MTF_cmass = intercept+MTF_count*slope_species_count2cmass

    return MTF_cmass
    
def convert_coords_to_lpjml_grid(lat,lon):
    """
    Take point coordinates and convert them to to the closest grid cell with a resolution of 0.5°x0.5°.
    
    Use to optain model or climate gridecell closest to field locations.
    """
    coords = [lat, lon]
    for i in range(len(coords)):
        if (abs(coords[i]) - abs(int(coords[i]))) < 0.5:
            #print(abs(coords[i]))
            if coords[i] < 0:
                coords[i] = int(coords[i]) - 0.25
            else:
                coords[i] = int(coords[i]) + 0.25
        else:
            if coords[i] < 0:
                coords[i] = int(coords[i]) - 0.75
            else:
                coords[i] = int(coords[i]) + 0.75
    return coords[0], coords[1]

def extract_climate(lats_in, lons_in, array):
    # Define lats and lons
    lats = np.arange(89.75, -89.75, -0.5)
    lons = np.arange(-179.75, 179.75, 0.5)

    # Convert the coordinates to the lpj grid
    coords = convert_coords_to_lpjml_grid(lats_in, lons_in) # [0] lats; lons [1]

    # Determine the lat index
    lat_idx = np.where(lats == coords[0])

    # Determine the lon index
    lon_idx = np.where(lons == coords[1])

    clim = array[lat_idx[0][0]][lon_idx[0][0]]

    return [clim, [lat_idx,lon_idx]]
    
    
def StemCarbon_whole(Species, DBH):
    path_ch = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/chojnacky_et_al_2013/'
    file = 'Chojnacky2014UpdatedSpecies_table5_parameters.csv'
    table5 = pd.read_csv(path_ch+file, encoding='utf-8-sig')

    path = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/TRY/211003_stem_carbon_content_TRY_407_16909/'
    file = 'TRY_mean_stem_carbon_content_species.csv'
    stem_c_content = pd.read_csv(path+file, encoding = "ISO-8859-1")
    try:
        row = table5[table5.Taxa.str.contains(Species.split(' ')[0])]
        
        # If there are multiple rows, Chojnacky distinguished by wood density
        if len(row) != 1:
            # Load Chojnacky Table 2,3 and 4 to check for how species was classified in Taxa
            file = 'Chojnacky2014UpdatedSpecies_wood_specific_gravity.csv'
            table_2_3_4 = pd.read_csv(path_ch+file, encoding='utf-8-sig')
            row_table_2_3_4 = table_2_3_4[table_2_3_4['Genus and species'].str.contains(Species)]
            
            ## Incase there are subspecies involved
            if (len(row_table_2_3_4) > 1):
                if (row_table_2_3_4.Taxa.values[0] == row_table_2_3_4.Taxa.values[1]):
                    row_table_2_3_4 = row_table_2_3_4.iloc[0,:].to_frame().T
            
            # Use taxa classification to read out the right row
            row = table5[table5.Taxa.str.contains(row_table_2_3_4.Taxa.item())]
        
        if (DBH > int(row['DBH max (cm)'].values)) | (DBH < int(row['DBH min (cm)'].values)):
            print('DBH outside defined range!')
            
        else:
            dry_mass = chojnacky_equation(row.beta_0, row.beta_1, DBH).item() # kg
            frac_c_content = stem_c_content[stem_c_content.SpeciesName.str.contains(Species)]['StdValue'] # fraction
            
            # If there is only one entry
            if len(frac_c_content) == 1:
                frac_c_content = frac_c_content.item()
                
            elif len(frac_c_content) > 1:
                frac_c_content = stem_c_content[stem_c_content.SpeciesName.str.contains(Species) &
                                    ~stem_c_content.SpeciesName.str.contains(' x ')]['StdValue'] # exclude hybrids
                frac_c_content = frac_c_content.item()
                
            else: # if there are no entries
                frac_c_content = 0.5 # Species not included in the data, assigning 0.5 c content
                print(Species+' not included in TRY stem carbon dataset! Using 50% carbon content.')
            
            c_content = dry_mass * frac_c_content # kg C
        return c_content
    except:
        print('Species / taxa not included!')
        
def chojnacky_equation(beta_0, beta_1, DBH):
        return np.exp(beta_0)*DBH**beta_1
        
        
def fracSnag_remain_at_t(t,frac_snag_remain):
    """
    Calculate k for an observation time t at which a fraction of snags remained.
    k is the decay rate in a single exponential model of the form y = e^{-kt}.

    t                -> time in years after mortality of snags occured.
    frac_snag_remain -> fraction of snags remaining.
    """
    F = 100 *(1-frac_snag_remain**(1/t)) # Equation 3
    k = np.log(-100/(F-100)) # Equation 2
    return k
    
def mod_single_exponential(x, k):
    return np.exp(-k*x)

def mod_linear(x,a,b):
    return (a+x*b)
    
def mod_linear_origin(x,b):
    return (x*b)

def mod_logistic(x,a,b,k):
    #y = a/(a+np.exp(k*(x-b))) # works for many
    y = 1/(a+np.exp(k*(x-b))) # works for Parish
#def mod_logistic(x,a,b):
#    y = 1/(a+np.exp(b*x))
    return y
    
def annualFall2k(F):
    """
    Equation 2
    """
    k = np.log(-(100/(F-100)))
    return k
    
    
def ref_max_dict(ref):
    """
    This defines the end point of the time series of the fraction of snags remaining at time t in years.
    """
    ref_max_dict = {'Zarnoch2013SnagForest':15, 'Parish2010SnagColumbia':200, 'Everett1999SnagUSA':200,
                    'Lyon1977AttritionMontana':50, 'Cline1980SnagOregon':100, 'Ganey2012Trends2007':125,
                    'Hogg2015FactorsCanada':70, 'Landram2002DemographyCalifornia':200, 'Holeksa2006ModelingPoland':200,
                    'Wilson2005DynamicsForest':100, 'Aakala2010CoarseDynamics':100, 'Lee1998DynamicsForests':50,
                    'Storaunet2002TimeForest':100, 'Mitchell1998FallOregon':20, 'Moorman1999SnagPiedmont':10,
                    'Russell2006SnagLogging':40, 'Chambers2014SnagForests':15, 'Bull1983LongevityWoodpeckers':30,
                    'Vanderwel2006AnTransitions':160, 'Moroni2006DisturbanceForests':34, 'Angers2010SnagSpecies':50,
                    'Dunn2015TemporalUSA':30, 'Angers2011TreeStudy':30, 'DeLong2008TemporalColumbia':200,
                    'Bond-Lamberty2008DecompositionChronosequence':100, 'Conner2005TheTexas':100,
                    'Ritchie2014EstablishmentForest':300, 'Onodera2015DoJapan':80, 'Huggard1999StaticSnags':200,
                    'Rhoades2020SnagfallUSA':60, 'Garber2005SnagMaine':100, 'Storaunet2004HowMethods':100,
                    'Vanderwel2006SnagForests':100, 'Molinas-Gonzalez2017FallMountain':30,
                    'Grayson2019PersistenceUSA':150, 'Dunn2012TemporalCascades':40, 'Woolley2019BeyondUSA':45,
                    'Yamasaki2006SnagHardwoods':30, 'Sinclair2004PersistenceAustralia':100, 'Chambers2005PonderosaArizona':100,
                    'Mäkinen2006PredictingFinland':70,
                   }
    max_x_interp = ref_max_dict[ref]
    return max_x_interp

def derive_MTF(ss, ref):
    
    """
    Calculate MTFcount from time series input.
    
    When DBH classes are provided calculate MTFcount for each DBH class and weigh to species level using the number of dead standing trees by DBH class.
    """
    # Print header
    pref = '#   '+ref+'   #'
    brdr = '#'+('-'*(len(pref)-2))+'#'
    print(brdr)
    print(pref)
    print(brdr+'\n')
    


    # Load data
    data = pd.read_csv('210910_MTF_input.csv', encoding='utf-8-sig')

    # Plotting
    cp = sns.color_palette("Set2", 13)
     
    # Initiliase variables
    step_size      = 0.001
    xl             = np.arange(0, ref_max_dict(ref), step_size) # first guess of the time dimension; used for plotting only
    MTF_method_idx = 1


    #--- Categorize data
    y_unit              = ss['y.unit'].iloc[0]
    DBH_mean_unique     = ss.DBH_mean.unique()
    DBH_classes_unique  = ss.DBH_class.unique()
    species_unique      = ss.species.unique()
    site_name_unique    = ss.site_name.unique()

    # Determine which MTFs need to be calculated
    if ~ss.cmass_whole.isna().all(): # if cmass_whole is filled
        category = ['count', 'cmass']
    else: # Count and Mass are the same category bc they only need to weighed by tree count
        if y_unit == 'count':
            category = ['count']
        elif y_unit == 'mass':
            category = ['c']

    print('No. sites: {0}\nNo. species: {1} \nNo. DBH classes: {2}\n'.format(
        len(site_name_unique), len(species_unique), len(DBH_classes_unique)))

    #--- Loop over categories –>   a:'count' & 'mass'   |   b: 'cmass'
    for cg in range(len(category)):
        print('\nCalculating MTF'+category[cg]+'!\n')
        
        #--- Loop over sites
        for s in range(len(site_name_unique)):

            # Distinguish between multiple sites & 1 site only
            if len(site_name_unique) == 1:
                df_sites                 = ss[ss.site_name == site_name_unique.item()]
                site                     = site_name_unique.item()
                site_species_unique      = df_sites.species.unique()
                site_species_list_unique = [[i,j] for i,j in zip(site_species_unique ,df_sites['species_scientific'].unique().tolist())]
                site_ntree_spe           = df_sites.groupby(['species', 'species_scientific', 'DBH_class']).first().species_ntree.values
                site_DBH_spe             = df_sites.groupby(['species', 'species_scientific', 'DBH_class']).first().DBH_mean.values
                site_species_percent     = site_ntree_spe / site_ntree_spe.sum() * 100
            else:
                df_sites                 = ss[ss.site_name == site_name_unique[s]]
                site_species_unique      = df_sites.species.unique()
                site_species_list_unique = [[i,j] for i,j in zip(site_species_unique ,df_sites['species_scientific'].unique().tolist())]
                site_ntree_spe           = df_sites.groupby(['species', 'species_scientific', 'DBH_class']).first().species_ntree.values
                site_DBH_spe             = df_sites.groupby(['species', 'species_scientific', 'DBH_class']).first().DBH_mean.values
                site_species_percent     = site_ntree_spe / site_ntree_spe.sum() * 100
                site = site_name_unique[s]
                
            #- Initialise plot depending on number of species
            if len(site_species_unique) > 4:
                if (len(site_species_unique) > 8) & (len(site_species_unique) <= 12):
                    divisor = 3
                elif len(site_species_unique) > 12:
                    divisor = 4
                else:
                    divisor = 2

                fig, ax = plt.subplots(divisor, 4, figsize=(20,15), sharey = True, sharex = True)
                row = 0
            else:
                fig, ax = plt.subplots(1,len(site_species_unique), figsize=(20,5), sharey = True, sharex = True)
            ax_idx = 0

            # Initialise arrays
            MTF_species        = np.zeros((len(site_species_unique),))+np.nan
            cmass_whole_spe    = np.zeros((len(site_species_unique),))+np.nan
            dbh_ntree_site_uni = np.zeros((len(site_species_unique),))+np.nan
            DBH_mean_spe       = np.zeros((len(site_species_unique),))+np.nan  # simple mean
            DBH_quadmean_spe   = np.zeros((len(site_species_unique),))+np.nan  # quadratic mean
            
            #--- Loop over species
            for spe in range(len(site_species_unique)):

                # Distinguish between multiple species& 1 species only
                if len(site_species_unique) == 1: # Single species
                    df_spe = df_sites[df_sites.species == site_species_unique.item()]
                    if category[cg] == 'cmass':
                        cmass_whole_dc = df_spe.cmass_whole.unique()
                    else:
                        cmass_whole_dc = np.array([np.nan]) # set to NaN
                    species             = site_species_unique.item()
                    species_list_unique = site_species_list_unique[spe]
                    species_percent     = site_species_percent[spe]
                    ntree_spe           = site_ntree_spe
                    spe_DBH_ntree_uni   = df_spe.groupby(['DBH_class','DBH_mean'],
                                                       sort=False)['DBH_ntree'].first().values
                    spe_DBH_percent     = spe_DBH_ntree_uni / np.sum(spe_DBH_ntree_uni) * 100 # relative contribution of DBH class to species
                else:# Multiple species
                    df_spe = df_sites[df_sites.species == site_species_unique[spe]]
                    if category[cg] == 'cmass':
                        cmass_whole_dc = df_spe.cmass_whole.unique()
                    else:
                        cmass_whole_dc = np.array([np.nan]) # set to NaN
                    species             = site_species_unique[spe]
                    species_list_unique = site_species_list_unique[spe]
                    species_percent     = site_species_percent[spe]
                    ntree_spe           = site_ntree_spe[spe]
                    spe_DBH_ntree_uni   = df_spe.groupby(['DBH_class','DBH_mean'],
                                                       sort=False)['DBH_ntree'].first().values
                    spe_DBH_percent     = spe_DBH_ntree_uni / np.sum(spe_DBH_ntree_uni) * 100
                
                print('\n',species)
                spe_DBH_mean_unique  = df_spe.DBH_mean.unique()
                spe_DBH_class_unique = df_spe.DBH_class.unique()
                

                # Initialise array to hold MTF results
                MTF_species_dc = np.zeros((len(spe_DBH_mean_unique),))+np.nan


                #--- Loop over DBH classes
                for dc in range(len(spe_DBH_mean_unique)):

                    # Distinguish between multiple DBH classes & 1 DBH class only
                    if len(DBH_mean_unique) == 1:
                        df_dbh_c = df_spe[df_spe.DBH_mean == spe_DBH_mean_unique.item()]
                        dbh = spe_DBH_mean_unique.item()
                    else:
                        df_dbh_c = df_spe[df_spe.DBH_mean == spe_DBH_mean_unique[dc]]
                        dbh = spe_DBH_mean_unique[dc]
                        
                    x = df_dbh_c.x.dropna().values
                    y = df_dbh_c.y.dropna().values
                    
                    #print(y[-1])
                    
                    donot_refs      = ['Everett1999SnagUSA', 'Chambers2005PonderosaArizona','Wilson2005DynamicsForest',]
                    donot_species   = ['Douglas fir',      # Everett1999SnagUSA
                                       #'Engelmann spruce', # Everett1999SnagUSA
                                       'Ponderosa pine',   # Chambers2005PonderosaArizona
                                       'Red maple',        # Wilson2005DynamicsForest
                                       ]
                    donot_dbh_class = ['23<41','41<64',               # Everett1999SnagUSA
                                       '10-24.9','25-50.9','51-88',   # Chambers2005PonderosaArizona
                                       '<10','10-20','20-30','30-50', # Wilson2005DynamicsForest
                                        ]

                    print('DBH class: ',spe_DBH_class_unique[dc],'cm; ', species)
                    
                    if (y[-1]>0.1): # if more than 10% are still standing
                        
                        if (ref in donot_refs) & (species in donot_species) & (spe_DBH_class_unique[dc] in donot_dbh_class):
                            # Do linear interpolation here for special cases where the line looks linear
                                                   
                            #- Interpolate data linearly
                            f = interp1d(x = x,
                                         y = y,
                                         kind='linear', fill_value='extrapolate')
                            y_intrp = f(xl) # extrapolate data
                            y_intrp[y_intrp < 0] = 0 # Set negative values to 0

                        else:
                            print('Logistic regression')
                            popt, pcov = opt.curve_fit(mod_logistic, x, y,
                                                       maxfev= 100000
                                                      )
                            y_intrp = mod_logistic(xl, *popt) # model data logistically
                    
                        
                        # Where does y reach 1%?
                        idx = np.where((y_intrp>0.009) & (y_intrp<=0.01))[0] # taking the first number where this condition is true
                        #print(len(idx))
                        if (len(idx) == 0) & ((ref == 'Everett1999SnagUSA')):# | (ref == 'Wilson2005DynamicsForest')): # exception only for everett!
                            xs = np.arange(0, 2000, step_size)
                            y_intrp = mod_logistic(xs, *popt)
                            idx = np.where((y_intrp>0.009) & (y_intrp<=0.01))[0]
                            y_intrp_mtf = y_intrp[:idx[0]]
                            xs = np.linspace(0, xs[idx[0]], idx[0]) # generate new x values
                            y_intrp = y_intrp[:len(xl)] # for plotting only
                        else:
                            #print(y_intrp,idx)
                            y_intrp_mtf = y_intrp[:idx[0]] # resize the y array until y >= 1%
                            xs = np.linspace(0, xl[idx[0]], idx[0]) # generate new x values
                            
                        #print(len(y_intrp_mtf), idx[0], y_intrp_mtf[-1], y_intrp[-1], len(xs))
                        
                        # Derive MTF
                        MTF = MTF_snags(xs, y_intrp_mtf, idx[0])[MTF_method_idx]
                        MTF_species_dc[dc] = MTF
                        
                        if len(DBH_mean_unique) > 1:
                            store_MTF_species_DBH_class_from_k(ref, MTF_species_dc[dc], MTF_method_idx, site, category[cg], species_list_unique, spe_DBH_ntree_uni[dc], spe_DBH_mean_unique[dc],
                                                               np.nan, spe_DBH_percent[dc], len(spe_DBH_class_unique))
                    else:
                        #- Interpolate data
                        #print(len(y),len(x),xl)
                        f = interp1d(x = x,
                                     y = y,
                                     kind='linear', fill_value='extrapolate')
                        y_intrp = f(xl)
                        y_intrp[y_intrp < 0] = 0 # Set negative values to 0
                        #print(y_intrp)
                        # Where does y reach 1%?
                        idx = np.where((y_intrp>0.009) & (y_intrp<=0.01))[0] # taking the first number where this condition is true
                        
                        #if len(idx) == 0:
                        #    idx = [len(xl)-1]
                        #print(idx, len(y_intrp), y_intrp[-1], xl[-1])
                        xs = np.linspace(0, xl[idx[0]], idx[0]) # generate new x values
                        y_intrp_mtf = y_intrp[:idx[0]] # resize the y array until y >= 1%
                        
                        #- Calculate MTF of record
                        #print('ALL MTF results: ',MTF_snags(xs, y_intrp_mtf, idx[0]))
                        MTF = MTF_snags(xs, y_intrp_mtf, idx[0])[MTF_method_idx]
                        MTF_species_dc[dc] = MTF
                        if len(DBH_mean_unique) > 1:
                            store_MTF_species_DBH_class_from_k(ref, MTF_species_dc[dc], MTF_method_idx, site, category[cg], species_list_unique, spe_DBH_ntree_uni[dc], spe_DBH_mean_unique[dc], #
                                                               np.nan, spe_DBH_percent[dc], len(spe_DBH_class_unique))

                    #- Store in input data
                    data.loc[(data.Reference == ref) &
                             (data.site_name == site) &
                             (data.species == species) &
                             (data.DBH_mean == dbh), 'MTF_record'] = MTF
                    data.loc[(data.Reference == ref) &
                             (data.site_name == site) &
                             (data.species == species) &
                             (data.DBH_mean == dbh), 'MTF_record_method'] = 'T2'
                    #data.to_csv('210910_MTF_input.csv', index=False, encoding='utf-8-sig')

                    #--- Plot
                    if len(site_species_unique) > 4:
                        ax[row][ax_idx].plot(xl, y_intrp, color=cp[dc], label=spe_DBH_class_unique[dc])
                        ax[row][ax_idx].plot(df_dbh_c.x, df_dbh_c.y, 'o', color=cp[dc], label='__None__')
                    elif len(site_species_unique) == 1:
                        ax.plot(xl, y_intrp, color=cp[dc], label=spe_DBH_class_unique[dc])
                        ax.plot(df_dbh_c.x, df_dbh_c.y, 'o', color=cp[dc], label='__None__')
                    else:
                        ax[ax_idx].plot(xl, y_intrp, color=cp[dc], label=spe_DBH_class_unique[dc])
                        ax[ax_idx].plot(df_dbh_c.x, df_dbh_c.y, 'o', color=cp[dc], label='__None__')


                #--- Weigh MTF to species level
                print('Weigh to species level',all(np.isnan(cmass_whole_dc)))
                if all(np.isnan(cmass_whole_dc)):
                    if ((len(spe_DBH_mean_unique) == 1) & (len(species_unique) == 1)):#
                        print('Transfer to site level only. No DBH classes!')# Only transfer MTFs to site level
                        MTF_species[spe]        = MTF_species_dc[dc]
                        #print('DBH mean1', spe_DBH_mean_unique)
                        if spe_DBH_mean_unique == '':
                            DBH_mean_spe[spe]   = np.nan
                        else:
                            DBH_mean_spe[spe]   = spe_DBH_mean_unique
                        DBH_quadmean_spe[spe]   = np.nan
                        print(MTF_species[spe])
                        store_MTF_species_from_k(ref, MTF_species[spe], MTF_method_idx, site, category[cg], species_list_unique, ntree_spe, DBH_mean_spe[spe], DBH_quadmean_spe[spe], species_percent, dc+1)
                    else:
                        #print(MTF_species_dc, spe_DBH_ntree_uni)
                        MTF_species[spe]        = np.sum(MTF_species_dc * spe_DBH_ntree_uni)/np.sum(spe_DBH_ntree_uni)
                        dbh_ntree_site_uni[spe] = np.sum(spe_DBH_ntree_uni)
                        
                        #print('DBH values: ', spe_DBH_mean_unique)
                        if any(spe_DBH_mean_unique == ''):
                            DBH_mean_spe[spe]       = np.nan
                            DBH_quadmean_spe[spe]   = np.nan
                        else:
                            DBH_mean_spe[spe]       = np.sum(spe_DBH_mean_unique * spe_DBH_ntree_uni) / np.sum(spe_DBH_ntree_uni) # mean weightes DBH
                            #print('DBH mean', np.mean(spe_DBH_mean_unique), np.sum(spe_DBH_mean_unique * spe_DBH_ntree_uni) / np.sum(spe_DBH_ntree_uni)) # mean weightes DBH)
                            DBH_quadmean_spe[spe]   = quadratic_mean_diam_DBH_class(spe_DBH_mean_unique, spe_DBH_ntree_uni)
                    print(MTF_species[spe])
                    store_MTF_species_from_k(ref, MTF_species[spe], MTF_method_idx, site, category[cg], species_list_unique, ntree_spe, DBH_mean_spe[spe], DBH_quadmean_spe[spe],species_percent, dc+1)
                        
                    
                else:# cmass whole only
                    if ((len(spe_DBH_mean_unique) == 1) & (len(species_unique) == 1)): # SHOULD this even exist?!
                        print('Cmass whole single')
                        MTF_species[spe]      = np.sum(MTF_species_dc[0] * cmass_whole_dc[0])/np.sum(cmass_whole_dc[0])
                        cmass_whole_spe[spe]  = np.sum(cmass_whole_dc[0])
                        DBH_mean_spe[spe]     = np.sum(spe_DBH_mean_unique[0] * spe_DBH_ntree_uni[0]) / np.sum(spe_DBH_ntree_uni[0]) # mean weightes DBH
                        DBH_quadmean_spe[spe] = quadratic_mean_diam_DBH_class(spe_DBH_mean_unique[0], spe_DBH_ntree_uni[0])
                        store_MTF_species_from_k(ref, MTF_species[spe], MTF_method_idx , site, category[cg], species_list_unique, ntree_spe, DBH_mean_spe[spe], DBH_quadmean_spe[spe], species_percent, dc+1)
                    else:
                        print('Cmass whole multiple','MTF_species_dc', MTF_species_dc)
                        MTF_species[spe]      = np.sum(MTF_species_dc * cmass_whole_dc)/np.sum(cmass_whole_dc)
                        cmass_whole_spe[spe]  = np.sum(cmass_whole_dc)
                        DBH_mean_spe[spe]     = np.sum(spe_DBH_mean_unique * spe_DBH_ntree_uni) / np.sum(spe_DBH_ntree_uni) # mean weightes DBH
                        #DBH_mean_spe[spe]     = np.mean(spe_DBH_mean_unique)
                        DBH_quadmean_spe[spe] = quadratic_mean_diam_DBH_class(spe_DBH_mean_unique, spe_DBH_ntree_uni) # quadratic mean diameter
                        store_MTF_species_from_k(ref, MTF_species[spe], MTF_method_idx, site, category[cg], species_list_unique, ntree_spe, DBH_mean_spe[spe], DBH_quadmean_spe[spe], species_percent, dc+1)

                # Plotdecor and increase indeces
                if len(site_species_unique) > 4:
                    if  ax_idx == 3:
                        if len(site_name_unique) > 1:
                            ax[row][ax_idx].set_title((ref+'\nSite:'+site+'\t Species: '+species).expandtabs())
                        else:
                            ax[row][ax_idx].set_title(ref+':\n'+species)
                        ax[row][ax_idx].set_ylabel('Frac remain [-]')
                        ax[row][ax_idx].set_xlabel('Time since death [years]')
                        ax[row][ax_idx].legend()
                        row += 1
                        ax_idx = 0
                    else:
                        if len(site_name_unique) > 1:
                            ax[row][ax_idx].set_title((ref+'\nSite:'+site+'\t Species: '+species).expandtabs())
                        else:
                            ax[row][ax_idx].set_title(ref+':\n'+species)
                        ax[row][ax_idx].set_ylabel('Frac remain [-]')
                        ax[row][ax_idx].set_xlabel('Time since death [years]')
                        ax[row][ax_idx].legend()
                        ax_idx += 1
                elif len(site_species_unique) == 1:
                    if len(site_name_unique) > 1:
                        ax.set_title((ref+'\nSite:'+site+'\t Species: '+species).expandtabs())
                    else:
                        ax.set_title(ref+':\n'+species)
                    ax.set_ylabel('Frac remain [-]')
                    ax.set_xlabel('Time since death [years]')
                    ax.legend()
                else:
                    if len(site_name_unique) > 1:
                        ax[ax_idx].set_title((ref+'\nSite:'+site+'\t Species: '+species).expandtabs())
                    else:
                        ax[ax_idx].set_title(ref+':\n'+species)
                    ax[ax_idx].set_title(ref+':\n'+species)
                    ax[ax_idx].set_ylabel('Frac remain [-]')
                    ax[ax_idx].set_xlabel('Time since death [years]')
                    ax[ax_idx].legend()
                    ax_idx += 1

            #--- Weigh MTF to site level
            if np.isnan(np.sum(cmass_whole_spe)): # weigh by number of species
                if np.isnan(np.sum(dbh_ntree_site_uni)): # if there were no ntrees reported
                    if len(MTF_species)> 1:
                        import sys
                        sys.exit('More than 1 species in MTF species. Need to calculate mean?')
                    else:
                        MTF_site          = MTF_species[0]
                        DBH_mean_site     = site_DBH_spe[0]
                        DBH_quadmean_site = np.nan
                else:
                    MTF_site          = np.sum(MTF_species * dbh_ntree_site_uni)/np.sum(dbh_ntree_site_uni)
                    if np.sum(site_DBH_spe == '') > 0:
                        DBH_mean_site     = np.nan
                        DBH_quadmean_site = np.nan
                    else:
                        DBH_mean_site     = np.sum(site_DBH_spe * site_ntree_spe) / np.sum(site_ntree_spe) # mean weighted DBH
                        #print('\n Mean Site DBH',DBH_mean_site,'\n',site_ntree_spe,site_DBH_spe)
                        DBH_quadmean_site = quadratic_mean_diam_DBH_class(site_DBH_spe, site_ntree_spe)
                print('MTF_site', MTF_site)
            else: # Weigh by cmass whole
                MTF_site          = np.sum(MTF_species * cmass_whole_spe)/np.sum(cmass_whole_spe)
                DBH_mean_site     = np.sum(site_DBH_spe * site_ntree_spe) / np.sum(site_ntree_spe) # mean weighted DBH
                #print('\n Mean Site DBH',DBH_mean_site,'\n',site_ntree_spe,site_DBH_spe)
                DBH_quadmean_site = quadratic_mean_diam_DBH_class(site_DBH_spe, site_ntree_spe)
                print('MTF_site', MTF_site)
            
            # Store result
            store_MTF_site(df_sites,MTF_site,MTF_method_idx,site,category[cg],site_species_unique,DBH_mean_site,DBH_quadmean_site)
            print(site)
              
            plt.subplots_adjust(hspace=0.5)
            plt.show()
        
    return


    
def store_MTF_species_from_k(reference,
                             MTF,
                             MTF_method,
                             site_name,
                             category,
                             species_list,
                             species_ntree,
                             species_mean_dbh,
                             species_quadmean_dbh,
                             species_prc,
                             n_dbh_class):
    
    #--- Load data
    
    if (type(MTF_method) == float) | (type(MTF_method) == int): # this is needed for data from the raw input file
        MTF_method_list = ['T1', 'T2', 'T3', 'T4', 'T5']
        MTF_method      = MTF_method_list[MTF_method]
    else:
        MTF_method      = MTF_method
    reference = reference
    
    file_wsg         = 'Chojnacky2014UpdatedSpecies_wood_specific_gravity.csv'
    path_wsg         = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/chojnacky_et_al_2013/'
    chojnacky_wsg    = pd.read_csv(path_wsg+file_wsg, encoding='utf-8-sig')
    path             = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/TRY/211003_stem_carbon_content_TRY_407_16909/'
    file             = 'TRY_mean_stem_carbon_content_species.csv'
    stem_c_content   = pd.read_csv(path+file, encoding = "ISO-8859-1")
    path_wd          = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/TRY/211003_wood_density_TRY_4_16909/'
    file_wd          = 'TRY_mean_wood_density_species.csv'
    wood_density     = pd.read_csv(path_wd+file_wd, encoding = "ISO-8859-1")
    path_wdurability = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/USA/Oberle2018WhenForests/'
    file_wdurability = 'Oberle2018_SI_Dataset1.csv'
    wood_durability  = pd.read_csv(path_wdurability+file_wdurability, encoding = "ISO-8859-1")
    path_n_content   = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/TRY/220302_stem_nitrogen_406_19485/'
    file_n_conten    = 'TRY_mean_stem_nitrogen_content_species.csv'
    stem_n_content   = pd.read_csv(path_n_content+file_n_conten, encoding = "ISO-8859-1")
    
            
    # Categorize data
    species            = species_list[0]
    species_scientific = species_list[1]
    spe_DBH_ntree_uni  = species_ntree

    # Here individual species information
    species_prc = species_prc
    wsp_species = chojnacky_wsg.loc[chojnacky_wsg['Genus and species'] == (species_scientific),
                                    'Wood specific gravity']
    wdurability_species = wood_durability.loc[(wood_durability.Binomial == species_scientific), 'composite']
    
    # Wood specific gravity
    if wsp_species.empty:
        wsp_species = np.nan
    else:
        wsp_species = wsp_species.item()
    
    # Stem carbon content
    frac_c_content = stem_c_content[stem_c_content.SpeciesName == species_scientific]['StdValue']
    frac_c_content_std = stem_c_content[stem_c_content.SpeciesName == species_scientific]['StdValue STD'] # Standard deviation
    if frac_c_content.empty:
        frac_c_content     = 0.5
        frac_c_content_std = np.nan
    else:
        frac_c_content     = frac_c_content.item()
        frac_c_content_std = frac_c_content_std.item()
    
    # Stem nitrogen content
    frac_n_content     = stem_n_content[stem_n_content.SpeciesName == species_scientific]['StdValue']
    frac_n_content_std = stem_n_content[stem_n_content.SpeciesName == species_scientific]['StdValue STD'] # Standard deviation
    if frac_n_content.empty:
        frac_n_content     = np.nan
        frac_n_content_std = np.nan
    else:
        frac_n_content     = frac_n_content.item()
        frac_n_content_std = frac_n_content_std.item()
    
    # Wood density
    wood_density_spe = wood_density[wood_density.SpeciesName == species_scientific]['StdValue']
    if wood_density_spe.empty:
        wood_density_spe = np.nan
    else:
        wood_density_spe = wood_density_spe.item()
    
    # Wood durability
    if wdurability_species.empty:
        wdurability_species = np.nan
    elif len(wdurability_species) > 1:
        wdurability_species = wdurability_species.mean().item()
    else:
        wdurability_species = wdurability_species.item()
    
    #--- Store result
    cols = ['MTF','MTF basis','MTF method','Species','Species scientific','Species % of total', 'DBH mean', 'Quadratic mean DBH',
            'No. DBH classes','Wood specific gravity','Wood density','Wood durability','Carbon content', 'Carbon content STD', 'Nitrogen content',
            'Nitrogen content STD', 'Site name','Reference', 'MTF method notes']

    data_count = [MTF, category, MTF_method, species, species_scientific, species_prc, species_mean_dbh,
                  species_quadmean_dbh, n_dbh_class, wsp_species, wood_density_spe, wdurability_species, frac_c_content, frac_c_content_std,
                  frac_n_content, frac_n_content_std, site_name, reference,  np.nan]
    df_MTF = pd.DataFrame({i:j for i,j in zip(cols, data_count)}, index= [0])

    if os.path.isfile('211030_MTF_species.csv') == True:
        MTF_loaded = pd.read_csv('211030_MTF_species.csv', encoding='utf-8-sig')
        MTF_loaded.loc[MTF_loaded['Site name'].isnull(), ['Site name']] = ''
        if (MTF_loaded.Reference == reference).any():
            print('removing the old ENTRY', category, species)
            MTF_loaded = MTF_loaded[~((MTF_loaded.Reference == reference) &
                                      (MTF_loaded['MTF basis'] == category) &
                                      (MTF_loaded['Species'] == species) &
                                      (MTF_loaded['Site name'] == site_name)
                                     )]
            MTF_store = pd.concat([MTF_loaded, df_MTF])
            MTF_store[cols].to_csv('211030_MTF_species.csv', index=False, na_rep='NA', encoding='utf-8-sig')
        else:
            print('entry does not exist')
            pd.concat([MTF_loaded, df_MTF]).to_csv('211030_MTF_species.csv', index=False, na_rep='NA', encoding='utf-8-sig')
    else:
        print('Creating file!')
        df_MTF.to_csv('211030_MTF_species.csv', index=False, na_rep='NA', encoding='utf-8-sig')
    return
    
    
def store_MTF_species_DBH_class_from_k(reference,
                                       MTF,
                                       MTF_method,
                                       site_name,
                                       category,
                                       species_list,
                                       species_ntree,
                                       species_mean_dbh,
                                       species_quadmean_dbh,
                                       species_prc, # contribution to total species trees
                                       n_dbh_class):

    #--- Load data

    if (type(MTF_method) == float) | (type(MTF_method) == int): # this is needed for data from the raw input file
        MTF_method_list = ['T1', 'T2', 'T3', 'T4', 'T5']
        MTF_method      = MTF_method_list[MTF_method]
    else:
        MTF_method      = MTF_method
    reference = reference

    file_wsg         = 'Chojnacky2014UpdatedSpecies_wood_specific_gravity.csv'
    path_wsg         = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/chojnacky_et_al_2013/'
    chojnacky_wsg    = pd.read_csv(path_wsg+file_wsg, encoding='utf-8-sig')
    path             = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/TRY/211003_stem_carbon_content_TRY_407_16909/'
    file             = 'TRY_mean_stem_carbon_content_species.csv'
    stem_c_content   = pd.read_csv(path+file, encoding = "ISO-8859-1")
    path_wd          = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/TRY/211003_wood_density_TRY_4_16909/'
    file_wd          = 'TRY_mean_wood_density_species.csv'
    wood_density     = pd.read_csv(path_wd+file_wd, encoding = "ISO-8859-1")
    path_wdurability = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/USA/Oberle2018WhenForests/'
    file_wdurability = 'Oberle2018_SI_Dataset1.csv'
    wood_durability  = pd.read_csv(path_wdurability+file_wdurability, encoding = "ISO-8859-1")
    path_n_content   = '/Users/antje/Documents/LUND/1912_modelling_SWD/2008_Dead_wood_treatment_of_models_review/data/snag_fall_rates/TRY/220302_stem_nitrogen_406_19485/'
    file_n_conten    = 'TRY_mean_stem_nitrogen_content_species.csv'
    stem_n_content   = pd.read_csv(path_n_content+file_n_conten, encoding = "ISO-8859-1")

            
    # Categorize data
    species            = species_list[0]
    species_scientific = species_list[1]
    spe_DBH_ntree_uni  = species_ntree

    # Here individual species information
    species_prc = species_prc
    wsp_species = chojnacky_wsg.loc[chojnacky_wsg['Genus and species'] == (species_scientific),
                                    'Wood specific gravity']
    wdurability_species = wood_durability.loc[(wood_durability.Binomial == species_scientific), 'composite']

    # Wood specific gravity
    if wsp_species.empty:
        wsp_species = np.nan
    else:
        wsp_species = wsp_species.item()
    
    # Stem carbon content
    frac_c_content = stem_c_content[stem_c_content.SpeciesName == species_scientific]['StdValue']
    frac_c_content_std = stem_c_content[stem_c_content.SpeciesName == species_scientific]['StdValue STD'] # Standard deviation
    if frac_c_content.empty:
        frac_c_content     = 0.5
        frac_c_content_std = np.nan
    else:
        frac_c_content     = frac_c_content.item()
        frac_c_content_std = frac_c_content_std.item()
    
    # Stem nitrogen content
    frac_n_content     = stem_n_content[stem_n_content.SpeciesName == species_scientific]['StdValue']
    frac_n_content_std = stem_n_content[stem_n_content.SpeciesName == species_scientific]['StdValue STD'] # Standard deviation
    if frac_n_content.empty:
        frac_n_content     = np.nan
        frac_n_content_std = np.nan
    else:
        frac_n_content     = frac_n_content.item()
        frac_n_content_std = frac_n_content_std.item()
    
    # Wood density
    wood_density_spe = wood_density[wood_density.SpeciesName == species_scientific]['StdValue']
    if wood_density_spe.empty:
        wood_density_spe = np.nan
    else:
        wood_density_spe = wood_density_spe.item()
    
    # Wood durability
    if wdurability_species.empty:
        wdurability_species = np.nan
    elif len(wdurability_species) > 1:
        wdurability_species = wdurability_species.mean().item()
    else:
        wdurability_species = wdurability_species.item()

    #--- Store result
    cols = ['MTF','MTF basis','MTF method','Species','Species scientific','Species % of total', 'DBH mean', 'Quadratic mean DBH',
            'No. DBH classes','Wood specific gravity','Wood density','Wood durability','Carbon content', 'Carbon content STD', 'Nitrogen content',
            'Nitrogen content STD', 'Site name','Reference','MTF method notes']

    data_count = [MTF, category, MTF_method, species, species_scientific, species_prc, species_mean_dbh,
                  species_quadmean_dbh, n_dbh_class, wsp_species,wood_density_spe, wdurability_species, frac_c_content, frac_c_content_std,
                  frac_n_content, frac_n_content_std, site_name, reference,  np.nan]
    df_MTF = pd.DataFrame({i:j for i,j in zip(cols, data_count)}, index= [0])

    if os.path.isfile('211213_MTF_species_DBH_class.csv') == True:
        MTF_loaded = pd.read_csv('211213_MTF_species_DBH_class.csv', encoding='utf-8-sig')
        MTF_loaded.loc[MTF_loaded['Site name'].isnull(), ['Site name']] = ''
        if (MTF_loaded.Reference == reference).any():
            print('removing the old ENTRY', category, species)
            MTF_loaded = MTF_loaded[~((MTF_loaded.Reference == reference) &
                                      (MTF_loaded['MTF basis'] == category) &
                                      (MTF_loaded['Species'] == species) &
                                      (MTF_loaded['Site name'] == site_name) &
                                      (MTF_loaded['DBH mean'] == species_mean_dbh)
                                     )]
            MTF_store = pd.concat([MTF_loaded, df_MTF])
            MTF_store[cols].to_csv('211213_MTF_species_DBH_class.csv', index=False, na_rep='NA', encoding='utf-8-sig')
        else:
            print('entry does not exist')
            pd.concat([MTF_loaded, df_MTF]).to_csv('211213_MTF_species_DBH_class.csv', index=False, na_rep='NA', encoding='utf-8-sig')
    else:
        print('Creating file!')
        df_MTF.to_csv('211213_MTF_species_DBH_class.csv', index=False, na_rep='NA', encoding='utf-8-sig')
    return
    
    
    
    
def store_MTF_site(ss,
                   MTF,
                   MTF_method_idx,
                   site,
                   category,
                   site_species_unique,
                   DBH_mean_site,
                   DBH_quadmean_site):

    #--- Load data
    MTF_method = ['T1', 'T2', 'T3', 'T4', 'T5']
    reference = ss.Reference.unique()[0]
    print(reference)
    sites = pd.read_csv('201217_MTF_sites.csv')
    
    sites.loc[sites['Site name'].isnull(), ['Site name']] = ''
    print(site)
    
    #display(sites[(sites.Reference == reference)])
    df_site = sites[(sites.Reference == reference) & (sites['Site name'] == site)]
    
    # Categorize data
    species_scientific_unique = ss.species_scientific.unique()
    dbh_ntree_unique = ss.DBH_ntree.unique()
    speciesDOM = df_site['Dominant species'].item()
    speciesDOM_prc = df_site['Dominant species % of total'].item()


    #--- Store result
    cols = ['MTF','MTF basis','MTF method','Dominant species','Species % of total', 'DBH mean',
            'Quadratic mean DBH', 'No. species contributing to site MTF','Site name', 'Reference',
            'MTF method notes']

    data_count = [MTF,category,MTF_method[MTF_method_idx],speciesDOM,speciesDOM_prc,
                  DBH_mean_site, DBH_quadmean_site, len(site_species_unique),
                  site, reference, np.nan]
    df_MTF = pd.DataFrame({i:j for i,j in zip(cols, data_count)}, index= [0])

    import os.path
    if os.path.isfile('210108_MTF.csv') == True:
        MTF_loaded = pd.read_csv('210108_MTF.csv')
        MTF_loaded.loc[MTF_loaded['Site name'].isnull(), ['Site name']] = ''
        if (MTF_loaded.Reference == reference).any():
            print('removing the old ENTRY')
            MTF_loaded = MTF_loaded[~((MTF_loaded.Reference == reference) &
                                      (MTF_loaded['MTF basis'] == category) &
                                      (MTF_loaded['Site name'] == site)
                                     )]
        else:
            print('entry does not exist')
    MTF_store = pd.concat([MTF_loaded, df_MTF])
    MTF_store[cols].to_csv('210108_MTF.csv', index=False, na_rep='NA', encoding='utf-8-sig')
    return

def store_MTF_site_from_k(MTF,
                          category,
                          MTF_method,
                          reference,
                          site_name,
                          speciesDOM,
                          speciesDOM_prc,
                          DBH_mean_site,
                          DBH_quadmean_site,
                          N_species_contributing_to_site_MTF,
                          MTF_method_notes):
    cols = [
        'MTF',
        'MTF basis',
        'MTF method',
        'Reference',
        'Site name',
        'Dominant species',
        'Species % of total',
        'DBH mean',
        'Quadratic mean DBH',
        'No. species contributing to site MTF',
        'MTF method notes',
        ]
    data = [MTF,
            category,
            MTF_method,
            reference,
            site_name,
            speciesDOM,
            speciesDOM_prc,
            DBH_mean_site,
            DBH_quadmean_site,
            N_species_contributing_to_site_MTF,
            MTF_method_notes]

    df_MTF = pd.DataFrame({i:j for i,j in zip(cols, data)}, index= [0])

    import os.path
    if os.path.isfile('210108_MTF.csv') == True:

        MTF_loaded = pd.read_csv('210108_MTF.csv')
        MTF_loaded.loc[MTF_loaded['Site name'].isnull(), ['Site name']] = '' # needed for identification

        if (MTF_loaded.Reference == reference).any():
            print('removing the old ENTRY')

            MTF_loaded = MTF_loaded[~((MTF_loaded.Reference == reference) &
                                      (MTF_loaded['MTF basis'] == category) &
                                      (MTF_loaded['Site name'] == site_name)
                                    )]
        else:
            print('Entry does not exist')
            
        MTF_store = pd.concat([MTF_loaded, df_MTF])
        MTF_store[cols].to_csv('210108_MTF.csv', index=False, na_rep='NA')
        
    else:
        df_MTF[cols].to_csv('210108_MTF.csv', index=False, na_rep='NA')
    return


def quadratic_mean_diam_DBH_class(DBH_class_mean, DBH_ntree):
    """
    'DBH_class_mean'   --> list of DBH class mean values
    'DBH_ntree'        --> number of trees / DBH class
    """
    DBH_ntree = np.round(DBH_ntree).astype('int') # first round to next integer, then assign datatype
    qmean_diam = np.sqrt(np.sum(np.repeat(DBH_class_mean, DBH_ntree)**2)/np.sum(DBH_ntree))
    return qmean_diam

def ODR_linear(x,y, plot_xbounds):
    """
    plot_xbounds   -> List of lower and upper x axis bounds
    
    """
    import scipy.odr as odr
    from scipy import stats

    def f(B, x):
        '''Linear function y = m*x + b'''
        # B is a vector of the parameters.
        # x is an array of the current x values.
        # x is in the same format as the x passed to Data or RealData.
        #
        # Return an array in the same format as y passed to Data or RealData.
        return B[0]*x + B[1]

    linear = odr.Model(f)

    mydata = odr.Data(x, y, wd=1./np.power(np.std(x),2), we=1./np.power(np.std(y),2))

    myodr = odr.ODR(mydata, linear, beta0=[0, 2.])

    myoutput = myodr.run()

    ss_res = myoutput.sum_square
    ss_tot = np.sum((y-np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('R2=',np.round(r_squared,2))

    N = len(y)

    # Degrees of freedom are n_samples - n_parameters
    df = N - 2  # equivalently, df = odr_out.iwork[10]
    beta_0 = 0  # test if slope is significantly different from zero
    t_stat = (myoutput.beta[0] - beta_0) / myoutput.sd_beta[0]  # t statistic for the slope parameter
    p_val = stats.t.sf(np.abs(t_stat), df) * 2
    print('Recovered equation: y={:3.2f}x + {:3.2f}, t={:3.2f}, p={:.2e}'.format(myoutput.beta[0],
                                                                             myoutput.beta[1],
                                                                             t_stat, p_val))
    beta = myoutput.beta
    X = np.linspace(plot_xbounds[0],plot_xbounds[1],100)
    Y = f(beta, X)

    return (X, Y, beta, df, r_squared, t_stat, p_val)

def read_chelsa30s_climate(clim_arr,lat,lon,var=np.nan):
    """
    Take point coordinates and read climate from chelsa grid cell of 30 arc seconds.
    lat      --> latitude, decimal degrees
    lon      --> longitude, decimal degrees
    var      --> 'tas' or 'wind' - vars like precip and temp seasonality are treated the same
    clim_arr --> stored as '.mat' file
    
    returns tas [°C] or precip [mm]
    """
    scale      = 0.1
    tas_offset = -273.15
    scale_wind = 0.001
    
    coords = [lat, lon]
    
    # Built latitude index
    lat_len   = 20880
    lat_start = -89.9960
    lat_end   = 83.9957
    step      = (abs(lat_start) + lat_end)/lat_len
    lats      = np.arange(lat_start, lat_end, step)
    
    # Built longitude index
    lon_len   = 43200
    lon_start = -179.996
    lon_end   = 170.9957
    step      = (abs(lon_start) + lon_end)/lon_len
    lons      = np.arange(lon_start, lon_end, step)
    
    # make combined list for looping
    index_lists = [lats, lons]
    
    
    # choose latitude and logitude index --> minimise distance to original coordinate
    chelsa_coords = np.zeros((2,), dtype='int')
    
    for i in range(len(coords)):
        coordinate = coords[i]
        full_index = index_lists[i]
        # Choosing the next larger and next smaller index value
        idx_larger  = full_index[full_index > coordinate][0]
        idx_smaller = full_index[full_index < coordinate][-1]
        idx = np.array([idx_larger, idx_smaller])
        
        # Check which one of them has the minimal distance to the original coordinate
        arr = [abs(coordinate-idx_larger), abs(coordinate-idx_smaller)]
        condition = (arr == min(arr))
        chelsa_coords[i] = np.where(full_index == idx[condition])[0][0]
    
    
    
    # Read value from climate array
    clim = clim_arr[chelsa_coords[0]][chelsa_coords[1]]
    
    # Conversion
    if var == 'temp':
        clim = clim * scale + tas_offset
    elif var == 'wind':
        clim = clim * scale_wind
    else:
        clim = clim * scale # can also be used for seasonality
        
    return clim
    
def read_soilTemp_lembrecht(clim_arr,lat,lon,var=np.nan):
    """
    Take point coordinates and read SoilTemperature from Lembrects2021Globaltemperature grid cell of 30 arc seconds.
    lat      --> latitude, decimal degrees
    lon      --> longitude, decimal degrees
    var      --> only monthly requires CHELSA scale of 0.1 - bioclim variables are correct - see plots
    clim_arr --> stored as '.mat' file
    """
    scale      = 0.1
    coords = [lat, lon]
    
    # Built latitude index
    lat_len   = 21121
    lat_start = -88.0042
    lat_end   = 87.9958
    step      = (abs(lat_start) + lat_end)/lat_len
    lats      = np.arange(lat_start, lat_end, step)
    
    # Built longitude index
    lon_len   = 43201
    lon_start = -179.9958
    lon_end   = 180.0042
    step      = (abs(lon_start) + lon_end)/lon_len
    lons      = np.arange(lon_start, lon_end, step)
    
    # make combined list for looping
    index_lists = [lats, lons]
    
    
    # choose latitude and logitude index --> minimise distance to original coordinate
    soilTemp_coords = np.zeros((2,), dtype='int')
    
    for i in range(len(coords)):
        coordinate = coords[i]
        full_index = index_lists[i]
        # Choosing the next larger and next smaller index value
        idx_larger  = full_index[full_index > coordinate][0]
        idx_smaller = full_index[full_index < coordinate][-1]
        idx = np.array([idx_larger, idx_smaller])
        
        # Check which one of them has the minimal distance to the original coordinate
        arr = [abs(coordinate-idx_larger), abs(coordinate-idx_smaller)]
        condition = (arr == min(arr))
        soilTemp_coords[i] = np.where(full_index == idx[condition])[0][0]
    
    
    
    # Read value from climate array
    clim = clim_arr[soilTemp_coords[0]][soilTemp_coords[1]]
    
    # Conversion
    if var == 'monthly':
        clim = clim * scale
    else:
        clim = clim
        
    return clim
    
def read_WorldClim_climate(clim_arr,lat,lon,var=np.nan):
    """
    Take point coordinates and read climate from WorldClim grid cell of 30 arc seconds or the gloabl aridity index.
    lat      --> latitude, decimal degrees
    lon      --> longitude, decimal degrees
    var      --> 'aridity' or nothing
    clim_arr --> stored as '.mat' file
    returns tas [°C] or precip [mm] or ai [-]
    """
    
    scalar_ai = 0.0001 # Global Aridity Index scalar
    
    coords = [lat, lon]
    
    # Built latitude index
    lat_len   = 21600
    lat_start = -89.9958
    lat_end   = 89.9958
    step      = (abs(lat_start) + lat_end)/lat_len
    lats      = np.arange(lat_start, lat_end, step)
    
    # Built longitude index
    lon_len   = 43200
    lon_start = -179.9958
    lon_end   = 179.9958
    step      = (abs(lon_start) + lon_end)/lon_len
    lons      = np.arange(lon_start, lon_end, step)
    
    # make combined list for looping
    index_lists = [lats, lons]
    
    
    # choose latitude and logitude index --> minimise distance to original coordinate
    chelsa_coords = np.zeros((2,), dtype='int')
    
    for i in range(len(coords)):
        coordinate = coords[i]
        full_index = index_lists[i]
        # Choosing the next larger and next smaller index value
        idx_larger  = full_index[full_index > coordinate][0]
        idx_smaller = full_index[full_index < coordinate][-1]
        idx = np.array([idx_larger, idx_smaller])
        
        # Check which one of them has the minimal distance to the original coordinate
        arr = [abs(coordinate-idx_larger), abs(coordinate-idx_smaller)]
        condition = (arr == min(arr))
        chelsa_coords[i] = np.where(full_index == idx[condition])[0][0]
    
    # Read value from climate array
    clim = clim_arr[chelsa_coords[0]][chelsa_coords[1]]
    
    
    if var == 'aridity':
        clim = clim * scalar_ai
    
    return clim


def store_single_loc(country,
                    region,
                    biome,
                    site_name,
                    y_coords,
                    x_coords,
                    min_height,
                    min_diam,
                    max_diam,
                    speciesDOM,
                    speciesDOM_prc,
                    mortality_cause,
                    management,
                    mean_age,
                    method_snag_field_measurement,
                    tree_share_to_snag,
                    n_plots,
                    plot_size,
                    ntree_Total,
                    tree_density,
                    Remeasurement_interval_avg,
                    survey_duration_total,
                    survey_start_yr,
                    survey_end_yr,
                    TSD_determination_method,
                    Model_type_fracRemain_snag,
                    Model_type_fracRemain_snag_sigCov,
                    Model_type_fracRemain_snag_nsigCov,
                    elevation,
                    mean_tC,
                    mean_precip_mm,
                    reference):
                     
    #——— SINGLE Location !!!!
    # Site data
    # Load excel file for storage
    import os.path
    
    if pd.isnull(site_name):
            site_name = ''
            
    cols = ['Country',
            'Region',
            'Site name',
            'Biome',
            'Y coords',
            'X coords',
            'Minimum height',
            'Minimum diameter',
            'Maximum diameter',
            'Dominant species',
            'Dominant species % of total',
            'Dominant_mortality cause',
            'Management',
            'Mean stand age',
            'Method_snagfall_measurement',
            'Fraction of live trees to snags',  # rest goes to SWD
            'N_plots',
            'Plot_area',  # hectare
            'ntree_Total',
            'tree_density', # individuals / ha
            'Remeasurement_interval_avg', # years
            'survey_duration_total', # years
            'survey_start_yr', # year
            'survey_end_yr', # year
            'TSD_determination_method',
            'Model_type_fracRemain_snag',
            'Model_type_fracRemain_snag_sigCov',
            'Model_type_fracRemain_snag_nsigCov',
            'Elevation',
            'Paper_T°C',
            'Paper_precip_mm',
            'Reference']
    # Define next row
    next_row = {
                'Country': country,
                'Region': region,
                'Biome': biome,
                'Site name': site_name,
                'Y coords': y_coords,
                'X coords': x_coords,
                'Minimum height': min_height,
                'Minimum diameter': min_diam,
                'Maximum diameter': max_diam,
                'Dominant species': speciesDOM,
                'Dominant species % of total':speciesDOM_prc,
                'Dominant_mortality cause': mortality_cause,
                'Management': management,
                'Mean stand age': mean_age,
                'Method_snagfall_measurement': method_snag_field_measurement,
                'Fraction of live trees to snags': tree_share_to_snag,  # rest goes to SWD
                'N_plots': n_plots,
                'Plot_area': plot_size,  # hectare
                'ntree_Total': ntree_Total,
                'tree_density': tree_density, # individuals / ha
                'Remeasurement_interval_avg': Remeasurement_interval_avg, # years
                'survey_duration_total': survey_duration_total, # years
                'survey_start_yr':survey_start_yr, # year
                'survey_end_yr':survey_end_yr, # year
                'TSD_determination_method': TSD_determination_method,
                'Model_type_fracRemain_snag': Model_type_fracRemain_snag,
                'Model_type_fracRemain_snag_sigCov': Model_type_fracRemain_snag_sigCov,
                'Model_type_fracRemain_snag_nsigCov': Model_type_fracRemain_snag_nsigCov, # not significant co-variats
                'Elevation': elevation,
                'Paper_T°C': mean_tC,
                'Paper_precip_mm': mean_precip_mm,
                'Reference': reference
               }
        
        # Load file if it exists
    if os.path.isfile('201217_MTF_sites.csv') == True:
        print('removing the old ENTRY')
        sites = pd.read_csv('201217_MTF_sites.csv')
        sites.loc[sites['Site name'].isna(), 'Site name'] = ''
            
        if (sites.Reference == reference).any():
            sites = sites[~((sites.Reference == reference) &
                           (sites['Site name'] == site_name)
                           )]
        else:
            print('Entry does not exist')
    
    # Store entry
        sites = sites.append(next_row,ignore_index=True)
        sites = sites[cols]
        sites.to_csv('201217_MTF_sites.csv',index=False, na_rep='NA')
    else:
        sites = pd.DataFrame(next_row, index=[0])
        sites = sites[cols]
        sites.to_csv('201217_MTF_sites.csv',index=False, na_rep='NA')
    

    return sites


def read_WorldClim_climate_10m(clim_arr,lat,lon,var=np.nan):
    """
    Take point coordinates and read climate from WorldClim grid cell of 10 arc minutes.
    lat --> latitude, decimal degrees
    lon --> longitude, decimal degrees
    var --> 'tas' or 'precip'
    
    returns tas [°C] or precip [mm]
    """
    
    coords = [lat, lon]
    
    # Built latitude index
    lat_len   = 1080
    lat_start = -89.9167
    lat_end   = 89.9167
    step      = (abs(lat_start) + lat_end)/lat_len
    lats      = np.arange(lat_start, lat_end, step)
    
    # Built longitude index
    lon_len   = 2160
    lon_start = -179.9167
    lon_end   = 179.9167
    step      = (abs(lon_start) + lon_end)/lon_len
    lons      = np.arange(lon_start, lon_end, step)
    
    # make combined list for looping
    index_lists = [lats, lons]
    
    
    # choose latitude and logitude index --> minimise distance to original coordinate
    chelsa_coords = np.zeros((2,), dtype='int')
    
    for i in range(len(coords)):
        coordinate = coords[i]
        full_index = index_lists[i]
        # Choosing the next larger and next smaller index value
        idx_larger  = full_index[full_index > coordinate][0]
        idx_smaller = full_index[full_index < coordinate][-1]
        idx = np.array([idx_larger, idx_smaller])
        
        # Check which one of them has the minimal distance to the original coordinate
        arr = [abs(coordinate-idx_larger), abs(coordinate-idx_smaller)]
        condition = (arr == min(arr))
        chelsa_coords[i] = np.where(full_index == idx[condition])[0][0]
    
    
    
    # Read value from climate array
    clim = clim_arr[chelsa_coords[0]][chelsa_coords[1]]
        
    return clim

def best_model(var_idx, var_list, var_list_in_flat, site):
    # Re-create regression string
    if site:
        columns_from_data = ['MTF']+[var_list_in_flat[j] for j in var_idx ]
    else:
        columns_from_data = ['MTF']+[var_list[j] for j in var_idx]
    
    best_model = ''
    for c in range(len(columns_from_data)-1):
            best_model = best_model+columns_from_data[c+1]+' + '
    best_model = columns_from_data[0]+' ~ '+best_model[:-3]
    return best_model


def step_wise_multiple_linear_reg_spe(df_MTF, alpha_sig, cycles, var_list_in, site=False, ignore_intercept_pval=False):
    
    # Copy df and change incompatible column names
    df_MTF_lm_full = df_MTF.copy()

    # Randomly shuffle rows of the dataframe in order to avoid serial autocorrelation
    df_MTF_lm_full = df_MTF_lm_full.sample(frac=1).reset_index(drop=True) # frac=1 mean that I am returning all rows in random order



    #--- Variables I want to test for
    var_list_in     = var_list_in
    alpha_sig       = alpha_sig
    random_var_list = True
    cycles          = cycles
    
    # Flatten input varlist
    var_list_in_flat = [s for sublist in [i.split(' + ') for i in var_list_in] for s in sublist]

    #- Initialise
    r_squared_cmass   = np.empty((len(var_list_in),cycles))
    AIC_cmass         = np.empty((len(var_list_in),cycles))
    RMSE              = np.empty((len(var_list_in),cycles))
    Durbin_Watson     = np.empty((len(var_list_in),cycles))
    White_test_r      = np.empty((len(var_list_in),cycles))
    reg_equations     = np.empty((len(var_list_in),cycles),dtype='object')
    reg_params        = np.empty((len(var_list_in),cycles),dtype='object')
    reg_pvals         = np.empty((len(var_list_in),cycles),dtype='object')
    reg_var_keep      = np.empty((len(var_list_in),cycles),dtype='object')
    var_idx           = [[] for i in range(cycles)] # list of lists for each cycle
    df_list           = []
    best_model_arr    = np.empty((cycles,),dtype='object')
    best_model_aic    = np.empty((cycles,),dtype='object')
    best_model_rmse   = np.empty((cycles,),dtype='object')
    best_model_dwt    = np.empty((cycles,),dtype='object') # Durbin Watson test
    best_model_wtr    = np.empty((cycles,),dtype='object') # White test result for homoscedasticity
    best_model_r2     = np.empty((cycles,),dtype='object')
    best_model_params = np.empty((cycles,),dtype='object')
    best_model_res    = np.empty((cycles,),dtype='object')


    for ci in range(cycles):
    
        # Shuffle variables randomly
        #print(random_var_list)
        if random_var_list:
            rng = np.random.default_rng()
            if site:
                var_list = list(rng.choice(var_list_in_flat[2:],len(var_list_in_flat)-2, replace=False))
                #print(var_list)
            else:
                var_list = rng.choice(var_list_in,len(var_list_in), replace=False)

        param_list  = []
        result_list = []
        
        if site:
            idx_len = len(var_list)+1
        else:
            idx_len = len(var_list)
            
        
        for i in range(idx_len):
            
            if site:
                if i == 0:
                    var_idx[ci].append(0)
                    var_idx[ci].append(1)
                    columns_from_data = ['MTF']+[var_list_in_flat[j] for j in var_idx[ci]]
                    #print(var_idx, columns_from_data)
                else:
                    var_idx[ci].append(i+1)
                    columns_from_data = ['MTF']+[var_list_in_flat[j] for j in var_idx[ci]]
                    
                
            else:
                var_idx[ci].append(i)
                columns_from_data = ['MTF']+[var_list[j] for j in var_idx[ci]]
                

            df_MTF_lm = df_MTF_lm_full[columns_from_data]

            # Regression string from variables
            reg_string = ''
            for c in range(len(columns_from_data)-1):
                reg_string = reg_string+columns_from_data[c+1]+' + '
            reg_string = columns_from_data[0]+' ~ '+reg_string[:-3]


            # Fitting linear model
            res_c = smf.ols(formula = reg_string, data = df_MTF_lm).fit()
            
            # Get coefficient specific response and prediction for RMSE
            y     = np.exp(df_MTF_lm.dropna(subset = columns_from_data).MTF.values) # generate MTF subset based on current coefficients & convert to years
            y_hat = np.exp(res_c.predict()) # convert to MTF years
                       
            
            # Store data
            r_squared_cmass[i,ci] = np.round(res_c.rsquared_adj,2)
            AIC_cmass[i,ci]       = np.round(res_c.aic)
            RMSE[i,ci]            = np.round(rmse_smf(y,y_hat),2)
            Durbin_Watson[i,ci]   = np.round(durbin_watson_test(res_c.resid),2) # Durbin Watson Test for serial autocorrelation
            try:
                White_test_r[i,ci]    = np.round(het_white(res_c.resid, res_c.model.exog)[1],3) # Conduct white test for homoscedasticity of the resisuals
            except:
                print('White test failed, assigning -1')
                White_test_r[i,ci] = -1
            reg_equations[i,ci]   = reg_string
            param_list.append(res_c.params)
            result_list.append(res_c)
            
            if ignore_intercept_pval:
                pvalues = res_c.pvalues[1:]
            else:
                pvalues = res_c.pvalues

            # Decide whether to discard or keep variable
            if i > 0:
                # Remove last cvar if R2 is ≤ R2[t-1] & AIC is ≥ than ANY previous AIC
                if ((pvalues > alpha_sig).sum() == 0):
                    if (any(r_squared_cmass[i,ci] <= r_squared_cmass[:i,ci]) &
                        any(AIC_cmass[i,ci] >= (AIC_cmass[:i,ci]-2))):
                        var_idx[ci] = var_idx[ci][:-1]                        # Removing the last variable
                        reg_var_keep[i][ci] = 'n'
                    else:
                        # Test the VIF of all variables
                        vif = fcy.Linear_Reg_Diagnostic(res_c).return_vif_table()
                        vif_idx = vif[vif['VIF Factor']>5].index.to_list()[1:] # only keep variables other than the Intercept
                        if vif_idx: # empty lists are considered False in python
                            var_idx[ci] = var_idx[ci][:-1]                        # Removing the last variable
                            reg_var_keep[i][ci] = 'n'
                        else:
                            reg_var_keep[i][ci] = 'y'
                else:
                    var_idx[ci] = var_idx[ci][:-1]
                    reg_var_keep[i][ci] = 'n'
            else: # For the variable just test whether the variable is significant and if r2 is larger than 0
                if site:
                        reg_var_keep[0][ci] = 'y'
                        reg_var_keep[1][ci] = 'y'
                
                else:
                    if ((pvalues > alpha_sig).sum() == 0):
                        if r_squared_cmass[i,ci] > 0:
                            reg_var_keep[i][ci] = 'y'
                            
                        else:
                            var_idx[ci] = var_idx[ci][:-1]
                            reg_var_keep[i][ci] = 'n'
                    else:
                        var_idx[ci] = var_idx[ci][:-1]
                        reg_var_keep[i][ci] = 'n'

        df_cols = ['r2','AIC','RMSE','DWT','WTR','var','keep']
        if site:
            df_data = [r_squared_cmass[:,ci], AIC_cmass[:,ci], RMSE[:,ci], Durbin_Watson[:,ci], White_test_r[:,ci], var_list_in[:1]+var_list, reg_var_keep[:,ci]]
        else:
            df_data = [r_squared_cmass[:,ci], AIC_cmass[:,ci], RMSE[:,ci], Durbin_Watson[:,ci], White_test_r[:,ci], var_list, reg_var_keep[:,ci]]
        #print(df_cols, df_data)
        df = pd.DataFrame({col:rows for col,rows in zip(df_cols,df_data)},
                          index=reg_equations[:,ci])
        df_list.append(df)

        #--- Store best model, R2 and AIC
        best_model_str       = best_model(var_idx[ci], var_list, var_list_in_flat,site)
        best_model_arr[ci]   = best_model_str
        if (best_model_str != 'MTF ~ '):
            bm_idx                = df.index.get_loc(best_model_str)
            best_model_r2[ci]     = df.iloc[bm_idx,0]
            best_model_aic[ci]    = df.iloc[bm_idx,1]
            best_model_rmse[ci]   = df.iloc[bm_idx,2]
            best_model_dwt[ci]    = df.iloc[bm_idx,3]
            best_model_wtr[ci]    = df.iloc[bm_idx,4]
            best_model_params[ci] = param_list[bm_idx]
            best_model_res[ci]    = result_list[bm_idx] # array of all best model results (parameters & coefficients)
        else:
            continue
        
        
        

    # All best models for each cycle
    result_cols = ['Best model', 'R2', 'AIC', 'RMSE','DWT','WTR']
    result_list = [best_model_arr, best_model_r2, best_model_aic, best_model_rmse, best_model_dwt,best_model_wtr]
    df_best_model = pd.DataFrame({col:rows for col,rows in zip(result_cols,result_list)},
                                 index=np.arange(0,cycles))
    df_best_model.loc[:, 'Model_idx'] = np.arange(0,cycles) # preserve index
    
    
    # Select the best models from all cycles
    df_best_model_first         = df_best_model.groupby('Best model').first()
    df_best_model_count         = df_best_model.groupby('Best model').count()['R2'].to_frame()
    df_best_model_count.columns = ['count']
    df_best_model_unique        = pd.concat([df_best_model_first,df_best_model_count],axis=1)
    df_best_model_unique        = df_best_model_unique.sort_values(['AIC','count']) # Table of the best model runs with R2, AIC and count
    
    #--- Many of these models are equivalent... So here summarising the models
    df_best_model_unique_reset = df_best_model_unique.reset_index()

    def clean_reg_covars(x):
        x = np.array(x, dtype='object')
        x = x[x != '+']
        x = x[x != '~']
        x = x[x != 'MTF']
        x = sorted(x)
        return x


    df_best_model_unique_reset.loc[:, 'Covars_ordered'] = df_best_model_unique_reset['Best model'].str.split().apply(clean_reg_covars)
    df_best_model_unique_reset.loc[:, 'Covars_n']       = df_best_model_unique_reset.loc[:, 'Covars_ordered'].apply(lambda x: len(x))
    df_best_model_unique_reset.loc[:, 'Covars_ordered'] = df_best_model_unique_reset.loc[:, 'Covars_ordered'].apply(lambda x: ' '.join(x))

    df_unique_models                      = df_best_model_unique_reset.groupby(['Covars_ordered']).sum()
    df_unique_models['R2']                = df_best_model_unique_reset.groupby(['Covars_ordered']).first()['R2']
    df_unique_models['AIC']               = df_best_model_unique_reset.groupby(['Covars_ordered']).first()['AIC']
    df_unique_models['RMSE']              = df_best_model_unique_reset.groupby(['Covars_ordered']).first()['RMSE']
    df_unique_models['DWT']               = df_best_model_unique_reset.groupby(['Covars_ordered']).first()['DWT']
    df_unique_models['WTR']               = df_best_model_unique_reset.groupby(['Covars_ordered']).first()['WTR']
    df_unique_models['Covars_n']          = df_best_model_unique_reset.groupby(['Covars_ordered']).first()['Covars_n']
    df_unique_models['Model_idx']         = df_best_model_unique_reset.groupby(['Covars_ordered']).first()['Model_idx']
    df_unique_models.loc[:, 'Best model'] = df_best_model_unique_reset.groupby(['Covars_ordered']).first()['Best model']
    df_unique_models                      = df_unique_models.sort_values('R2', ascending = False)
    
    
    display(df_unique_models)
    print(df_unique_models.to_latex()) # print out the result
    
    
    #-----------------#
    #   Provisional   #
    #-----------------#
    # Selecting the best model
    # --> the least amount of variables and the highest number of occurrences for the highest R2
    df_best_model_unique_t = df_unique_models.reset_index()
    df_best_model_unique_t = df_best_model_unique_t[(df_best_model_unique_t.R2 == df_best_model_unique_t.R2.max())]
                                                      
                                                      
    THE_best_model = df_best_model_unique_t[df_best_model_unique_t['count'] == df_best_model_unique_t['count'].max()]
        
    display(THE_best_model)

    bm_df_idx = THE_best_model['Model_idx'].item()

    res_export = best_model_res[bm_df_idx] # the result object of the final best model
    
    return df_unique_models, best_model_res, res_export
    
    
def read_SoilWatergrids(clim_arr,lat,lon):
    """
    Take point coordinates and read climate from SoilWatergridsv1 grid cell of 30 arc seconds.
    lat      --> latitude, decimal degrees
    lon      --> longitude, decimal degrees
    clim_arr --> soil water grids mean annual soil water saturation of the frist layer, stored as .mat file
    """
      
    coords = [lat, lon]
    
    # Built latitude index
    lat_len   = 720
    lat_start = -89.875
    lat_end   = 89.875
    step      = (abs(lat_start) + lat_end)/lat_len
    lats      = np.arange(lat_start, lat_end, step)
    
    # Built longitude index
    lon_len   = 1440
    lon_start = -179.875
    lon_end   = 179.875
    step      = (abs(lon_start) + lon_end)/lon_len
    lons      = np.arange(lon_start, lon_end, step)
    
    # make combined list for looping
    index_lists = [lats, lons]
    
    
    # choose latitude and logitude index --> minimise distance to original coordinate
    swg_coords = np.zeros((2,), dtype='int')
    
    for i in range(len(coords)):
        coordinate = coords[i]
        full_index = index_lists[i]
        # Choosing the next larger and next smaller index value
        idx_larger  = full_index[full_index > coordinate][0]
        idx_smaller = full_index[full_index < coordinate][-1]
        idx = np.array([idx_larger, idx_smaller])
        
        # Check which one of them has the minimal distance to the original coordinate
        arr = [abs(coordinate-idx_larger), abs(coordinate-idx_smaller)]
        condition = (arr == min(arr))
        swg_coords[i] = np.where(full_index == idx[condition])[0][0]
    
    # Read value from climate array
    clim = clim_arr[swg_coords[0]][swg_coords[1]]
    
    return clim
    
    
def climate_decomposition_index(T, PPT, PET):
    """
    Input values need to be monthly!
    From Adair2008SimpleClimates
    
    T   --> mean monthly airtemperature
    PET --> potential evaportranspiration (ET0)
    PPT --> precipitation
    
    Lloyd and Taylor (1994) variable Q10 temperature function (Ft).
    Source: https://adairlab.weebly.com/how-to-calculate-cdi.html
    """
    
    Ft = 0.5766 * np.exp(308.56 * ((1/56.02)-(1/((273+T)-227.13))))
    
    Fw = 1/ (1+30 * np.exp(-8.5 * (PPT/PET)));

    cdi = Ft * Fw

    return cdi
    
    
    
def replace_conversion_slope(path_filename, find_string, add_number2string):
    """
    Open a file and replace a string.
    Here used to update conversion slope between MTFcount and MTF size.
    """
    
    new_line = find_string + str(add_number2string)+'            \n'
    
    with open(path_filename, 'r') as file :
        filedata = file.read()

        if filedata.find(find_string) == -1:
            raise ValueError('Search string not in file!')

    with open(path_filename, 'r') as myFile:
        
        
        for num, line in enumerate(myFile, 1):
            if find_string in line:
                print('found at line:', num)
                print(line)
                pattern = line

    
    # Set in new line at the old line (pattern)
    f = open(path_filename,'r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace(pattern,new_line)

    f = open(path_filename,'w')
    f.write(newdata)
    f.close()
    
    return None
    
    
def classify_mortality(mtf):
    mtf_mort               = mtf.dropna(subset=['Dominant_mortality cause']).copy()
    mtf_fire_killed        = mtf_mort[mtf_mort['Dominant_mortality cause'].str.contains('fire') |
                               mtf_mort['Dominant_mortality cause'].str.contains('Fire')
                              ].copy()
    ref_mtf_fire_killed    = list(mtf_fire_killed.Reference.unique())
    mtf_beetle_killed      = mtf_mort[mtf_mort['Dominant_mortality cause'].str.contains('beetle') |
                                 mtf_mort['Dominant_mortality cause'].str.contains('moth')
                                ].copy()
    ref_mtf_beetle_killed  = list(mtf_beetle_killed.Reference.unique())
    mtf_drought_killed     = mtf_mort[mtf_mort['Dominant_mortality cause'].str.contains('Drought')].copy()

    ref_fire_beetle_killed = ref_mtf_fire_killed+ref_mtf_beetle_killed

    print('{} references reported fire as dominant mortality cause.'.format(len(ref_mtf_fire_killed)))
    print('{} references reported beetles as dominant mortality cause.'.format(len(ref_mtf_beetle_killed)))
    print('{} references reported drought as dominant mortality cause.'.format(len(mtf_drought_killed.Reference.unique())))
    print(len(ref_fire_beetle_killed))
    mtf.loc[:, 'Beetle']                        = mtf.Reference.apply(lambda x: 1 if x in ref_mtf_beetle_killed else 0)
    mtf.loc[:, 'Mort_other']                    = mtf.Reference.apply(lambda x: 1 if x not in ref_fire_beetle_killed else 0)
    mtf.loc[(mtf.Fire == 1), 'Mortality']       = 'Fire'
    mtf.loc[(mtf.Beetle == 1), 'Mortality']     = 'Insects'
    mtf.loc[(mtf.Mort_other == 1), 'Mortality'] = 'Other'
    mtf.loc[(mtf.Mortality.isna()),'Mortality'] = 'Other'
    return mtf
    
    
def regression_wrapper(df, model_list, alpha, mtf_type, climate_data, standardise_covars):
    
    # Number of models
    n_models = len(model_list)
    
    # Copy df and change incompatible column names
    df_MTF_lm = df.copy()
    
    #- Initialise
    r_squared        = np.empty((n_models))
    AIC              = np.empty((n_models))
    RMSE             = np.empty((n_models))
    Durbin_Watson    = np.empty((n_models))
    White_test_r     = np.empty((n_models),dtype='object')
    reg_params       = np.empty((n_models),dtype='object')
    reg_VIF_ok       = np.empty((n_models),dtype='object') # Is VIF of any regression ≥ 5?
    condition_number = np.empty((n_models),dtype='object') # Condition number valid?
    eigenvals        = np.empty((n_models),dtype='object') # Eigenvalues valid?
    pvals_pass       = np.empty((n_models),dtype='object')
    pvals_notpass_p  = np.empty((n_models),dtype='object')
    n_obs            = np.empty((n_models))
    param_list       = []
    result_arr       = np.empty((n_models),dtype='object')

    # Loop over models
    for m in range(n_models):
        
        # Fitting linear model
        res_c = smf.ols(formula = model_list[m], data = df_MTF_lm).fit()
        #print(model_list[m])
        
        # Get coefficient specific response and prediction for RMSE
        y     = np.exp(df_MTF_lm.MTF.values) # generate MTF subset based on current coefficients & convert to years
        y_hat = np.exp(res_c.predict())      # convert to MTF years
        
        
        
        #--- Store regression metrics
        r_squared[m]        = np.round(res_c.rsquared_adj,2)
        AIC[m]              = np.round(res_c.aic,1)
        RMSE[m]             = np.round(rmse_smf(y,y_hat),2)
        Durbin_Watson[m]    = np.round(durbin_watson_test(res_c.resid),2) # Durbin Watson Test for serial autocorrelation
        try:
            White_test_r[m] = (het_white(res_c.resid, res_c.model.exog)[1] < 0.05) == False  # Conduct white test for homoscedasticity of the resisuals
        except:
            print('White test failed, assigning -1')
            White_test_r[m] = -1
        condition_number[m] = res_c.condition_number <= 200          # Needs to be True 1000
        eigenvals[m]        = (res_c.eigenvals[-1] < 1e-10) == False # Needs to be True for valid model
        n_obs[m]            = res_c.nobs
        pvals_pass[m]       = any((res_c.pvalues[1:] > alpha).values) == False
        
        # Insignificant params
        df = (res_c.pvalues[1:] > alpha).to_frame().reset_index()
        pvals_notpass_p[m]  = ', '.join(df.loc[df[0] == True, 'index'].tolist())
        
        param_list.append(res_c.params)
        result_arr[m]       = res_c
        
        # Test the VIF of all variables
        try:
            vif = fcy.Linear_Reg_Diagnostic(res_c).return_vif_table()
            #display(vif)
            vif_idx = vif[vif['VIF Factor']>=5].index.to_list()[1:] # only consider variables other than the Intercept
            if vif_idx: # empty lists are considered False in python
                reg_VIF_ok[m] = 'fail'
            else:
                reg_VIF_ok[m] = 'pass'
        except:
            reg_VIF_ok[m] = 'fail'
        
        # Result table
        df_cols = ['r2','AIC', 'RMSE','DWT','WTR','VIF','CondN','Eigvals','prms_sig','prms_fail','N_obs', 'ID','climate','standardised_covars']
        df_data = [r_squared, AIC, RMSE, Durbin_Watson, White_test_r,
                   reg_VIF_ok,condition_number,eigenvals,pvals_pass,pvals_notpass_p,n_obs,np.arange(0,n_models), climate_data, standardise_covars]
        df_results = pd.DataFrame({col:rows for col,rows in zip(df_cols,df_data)},
                          index=model_list)
        df_results['MTF_basis'] = mtf_type
        
        # Model selection
        df_results = df_results[(df_results.WTR != False)]
        df_results = df_results[(df_results.VIF != 'fail')]
        df_results = df_results[(df_results.CondN != False)]
        df_results = df_results[(df_results.Eigvals != False)]
        
        # Compute ∆ AIC as ∆i = AICi - AIC_min
        df_results['d_AIC'] = abs(df_results.AIC.min() - df_results.AIC)
        df_result_ss = df_results.loc[:,['r2','AIC','d_AIC','RMSE','DWT','prms_sig','prms_fail','N_obs','ID','MTF_basis','climate','standardised_covars']].copy()
        df_result_ss = df_result_ss.sort_values(['r2','AIC','RMSE'],ascending=[False,True,True]).copy()
        
        # Get subset of the result objects
        result_arr_ss = result_arr[df_result_ss.ID.tolist()]
        
        # Reset the ID in the table to match result object index
        df_result_ss.ID = np.arange(0,len(df_result_ss))
        
    return df_result_ss, result_arr_ss
    
    
    
def MTF_database_import(file, sheet_name, engine=np.nan):
    """
    Import MTF database and set MultiIndex
    """
    if pd.isna(engine):
        df_raw = pd.read_excel(file, sheet_name=sheet_name, header=[0,1], )
    else:
        df_raw = pd.read_excel(file, sheet_name=sheet_name, header=[0,1], engine='openpyxl')
    
    cols = df_raw.columns.to_list()
    cols = [list(i) for i in (cols)]

    for i,i_list in enumerate(cols):
        for j,j_item in enumerate(i_list):
            if 'Unnamed' in j_item:
                cols[i][j] = ''

    cols = np.array([np.array(i) for i in cols])
    df_raw.columns = pd.MultiIndex.from_arrays(cols.T)
    if pd.isna(engine):
        df = df_raw.iloc[1:,1:].reset_index(drop=True)
    else:
        df = df_raw.iloc[:,1:].reset_index(drop=True)
    return df
