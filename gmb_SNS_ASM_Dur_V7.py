# Auxiliary code for __
# This code generats results figures for the study mentioned above and some
# supplementary material

# This code has been developed on Python 3.12.5 with Spyder 6
# The code has been designed to be abel to run by sections on Spyder, other
# interpreters may not identify the same section berakpoints.

# Further details on readme document

#%% Initialize workspace
# Import required libraries
import os
import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt
import seaborn           as sns

from scipy    import stats

import statsmodels.formula.api as smf

# Figure generation workspace
# Uncomment the following 6 lines to obtain better quality figure generation
import matplotlib.font_manager as ft
import matplotlib
matplotlib.use('cairo')
plt.rc('font', family='arial')
plt. rcParams.update({'font.size':9})
font = ft.FontProperties(family = 'arial')

# Define the Origin path for all data
path_Org = os.getcwd()

# Clustering data importins
sz_Clust = pd.read_csv(path_Org + '/data/rednmf_l10.csv')
Stt_Freq_asm = pd.read_csv(path_Org + '/data/rednmf_l10_szFrq.csv')

# Figure saving flag: True = Save figures
save_fig = True

# Figure saving output folders
path_save_fig = path_Org + '/ASM_Result_figs/'

# Number of decimas on statistical analysis show
Stt_Dgt = 4

# SVG/PDF figure folders
if not os.path.exists(path_save_fig + '/PDF'):
    os.makedirs(path_save_fig + '/PDF')

#%% Results: Seizure dutation and Frecuancy with ASM levels

# Duration for all seizures available on the dataset
rsp = "Duration_10"

# Figure Aperance variables
Y_Lim1 = [0,3]
Y_Lim2 = [0,3]
Y_Label= "Duration (log10)"

X_Lim_asm = [-10,1]
X_Lim_tod = [0,24]

palt = "mako"
Hue_Patt = ["Normo","Under"]
Axs_Patt = ["Under","Normo"]
Cater = "ASM_stt"

# MELM formula for duation against ASM levels and TofD
# sTOD+cTOD     -> Incorporation ciradian rythm (24h) as a circular variable
# sTOD_h+cTOD_h -> Incorporation ultradian 12h rythm as a circular variable
prd = "ASM_lvl + sTOD+cTOD + sTOD_h+cTOD_h"

# Clear dataset from non-relevant data
rm_dx = np.logical_or(np.isnan(sz_Clust[rsp]),np.isinf(sz_Clust[rsp]))
stData = sz_Clust.drop(np.where(rm_dx)[0]).reset_index(drop=True)

# Fit the mixed effect model with patietn identifier group-level effect
mdl_Chr_Ams = smf.mixedlm(rsp + " ~ " + prd , stData, groups=stData.patient_id).fit()
mdl_Chr_Ams.summary()

# Relevent effect sizes and p-values for figure generation
rnd_eff = mdl_Chr_Ams.random_effects
int_eff = mdl_Chr_Ams.params.Intercept

asm_eff = mdl_Chr_Ams.params.ASM_lvl
asm_pvl = mdl_Chr_Ams.pvalues.ASM_lvl

plt.figure(figsize=[15,5])

# Plot raw data for 3 individuals: ASM level vs Seizure duration --------------
# Patient Ids to show
id_rep = [78,46,48]

# Color code
id_col = [[59/255, 117/255, 175/255],
          [235/255,125/255, 48/255],
          [80/255, 156/255, 61/255]]

# Plot data poitns and Fixed effect individually
for i in range(0,len(id_rep)):
    plt.subplot(1,3,2)
   
    # Identify subject data
    idx = np.where(sz_Clust.patient_id == id_rep[i])[0]
    
    # Data poitns scatterplot
    sns.scatterplot(x=sz_Clust.ASM_lvl[idx], y=sz_Clust[rsp][idx], 
                    label="id" + str(id_rep[i]), color = id_col[i])
    
    # Generate regresion with MELM data
    grp_eff = rnd_eff.get(id_rep[i]).Group     # Group effect
    X = np.linspace(min(sz_Clust.ASM_lvl[idx]),# ASM level range 
                    max(sz_Clust.ASM_lvl[idx]))
    Y = X * asm_eff + grp_eff + int_eff        # Regresion
    sns.lineplot(x=X, y=Y, color = id_col[i])  # Line with matching color code

# Mark the ASM levels correponding to NOrmal-dose in grey
plt.fill_between(range(X_Lim_asm[0],X_Lim_asm[1]+1), Y_Lim1[0], Y_Lim1[1], 
                 np.array(range(X_Lim_asm[0],X_Lim_asm[1]+1))>=-1, 
                 alpha=0.1, color='k')    

plt.subplot(1,3,1)
plt.xlim(X_Lim_asm)
plt.ylim(Y_Lim1)

plt.legend()
plt.ylabel(Y_Label)
plt.xlabel("normalized ASM plasma concentration")

plt.title("Three Example Patients")

# Plot all individuals: ASM level vs Seizure duration with MELM ---------------
plt.subplot(1,3,3)

# Remove the group effect to all subjects
FL = np.repeat(np.nan, stData.shape[0])
for i in rnd_eff.keys():
    dt_dx = stData.patient_id == i
    if any(dt_dx):
        FL[dt_dx] = stData[rsp][dt_dx] - rnd_eff.get(i).Group
stData[rsp] = FL

# Plot the datapoints color codding the ASM_level [Normal-dose | Under-dose] 
sns.scatterplot(data = stData, x="ASM_lvl", y=rsp, 
                hue = Cater, palette = palt, hue_order = Hue_Patt)

# Generate regresion for all datapoints
X = np.linspace(X_Lim_asm[0], X_Lim_asm[1])
Y = int_eff
Y += X * asm_eff * (asm_pvl < 0.05)
sns.lineplot(x=X, y=Y, color = "red", label = "ASM lvl. - Fixed eff.")

# Mark the normal-dose period
plt.fill_between(range(X_Lim_asm[0],X_Lim_asm[1]+1), Y_Lim2[0], Y_Lim2[1], 
                 np.array(range(X_Lim_asm[0],X_Lim_asm[1]+1))>=-1, 
                 alpha=0.1, color='k')

plt.xlim(X_Lim_asm)
plt.ylim(Y_Lim2)
plt.ylabel(Y_Label)
plt.xlabel("normalized ASM plasma concentration")

# Include sample size and ASM level fixed effect statistics
plt.title("All Seizures (" +
          str(stData.shape[0])+ " Szs. | "+
          str(len(np.unique(stData.patient_id)))+ " Pts.)" +
          "\n MELM - ASM lvl. eff=" + str(round(asm_eff,Stt_Dgt)) + 
          ", p=" + str(round(asm_pvl,Stt_Dgt)))

# Plot seizure frecuancy changes between ASM level ranges on all subjects
axs = plt.subplot(1,3,1)

col_thr = 0
for i in range(0,Stt_Freq_asm.shape[0]):
    # Colorcode the data based on difference in seizure frecuancy
    if np.diff(Stt_Freq_asm.values[i])>col_thr: # Increase -> in blue
        COL = 'blue'
        STL = ':'
        MRK = 'o'
        ALP = 1
    elif np.diff(Stt_Freq_asm.values[i])<-col_thr: # Decrease -> in red
        COL = 'red'
        STL = '-.'
        MRK = '^'
        ALP = 1
    else: # No Change -> in black
        COL = 'black'
        STL = '--'
        MRK = 'X'
        ALP = 1
        
    # Plot a line per subject 
    # Square root of data for better visualization with Y^2 axis
    sns.lineplot(y=np.sqrt(Stt_Freq_asm.values[i][::-1]),x=["Under","Normo"],
                 color=COL, alpha = ALP, legend = False, linestyle = STL,
                 marker = MRK)

xl = np.array([-0.5,0,0.5,1,1.5])
yl = [0,7.5]

# Indicate the Normal-dose data with grey background
plt.fill_between(xl, yl[0], yl[1], xl>=0.5, alpha=0.1, color='k')
plt.xlim([xl[0],xl[-1]])

# Square Axis
y_tick = np.linspace(yl[0],round(yl[1]),round(yl[1])+1)
plt.yticks(y_tick, np.square(y_tick).astype(str))

plt.ylabel("sz/day")
plt.ylim(yl)

# Statistical comparaison Comparaison
Freq_SignRank = stats.wilcoxon(Stt_Freq_asm.Under,Stt_Freq_asm.Normo, method = 'asymptotic', 
                               alternative = 'greater')

R_stat = Freq_SignRank.zstatistic/np.sqrt(Stt_Freq_asm.shape[0])

# Include sample size and statistics
plt.title("Sz. Frequency[" + str(Stt_Freq_asm.shape[0]) + 
          "] \nSinged-Rank R=" + str(round(R_stat,Stt_Dgt)) + 
          ", p="+ str(round(Freq_SignRank.pvalue,Stt_Dgt)))

if save_fig:
    plt.savefig(path_save_fig + "/RLT_fDS_Durt_Freq.svg", format='svg', 
                dpi=3000, bbox_inches='tight')
    plt.savefig(path_save_fig + "/PDF/RLT_fDS_Durt_Freq.pdf", format='pdf', 
                dpi=3000, bbox_inches='tight')

#%% Results: Seizure dutation, Seizure Staes and ASM levels

# Full log10 duration for seizures within the Seizure Sates analysis
rsp = "DSS_Full"

# Figure Aperance variables
Y_Lim1 = [0,3]
Y_Lim2 = [0,3]
Y_Label= "Duration (log10)"

X_Lim_asm = [-10,1]
X_Lim_tod = [0,24]

palt = "mako"
Hue_Patt = ["CMN_Normo","CMN_Under","UND_Under"]
Axs_Patt = ["UND_Under","CMN_Under","CMN_Normo"]
Cater = "Sz_Class"

# MELM for for duation against ASM levels and TofD
prd = "ASM_lvl + sTOD+cTOD + sTOD_h+cTOD_h"

# Clear dataset from non-relevant data
rm_dx = np.logical_or(np.isnan(sz_Clust[rsp]),np.isinf(sz_Clust[rsp]))
stData = sz_Clust.drop(np.where(rm_dx)[0]).reset_index(drop=True)

# Fit the mixed effect model with patietn identifier group-level effect
mdl_Chr_SSS = smf.mixedlm(rsp + " ~ " + prd , stData, groups=stData.patient_id).fit()
mdl_Chr_SSS.summary()

# Relevent effect sizes and p-values for figure generation
rnd_eff = mdl_Chr_SSS.random_effects
int_eff = mdl_Chr_SSS.params.Intercept

asm_eff = mdl_Chr_SSS.params.ASM_lvl
asm_pvl = mdl_Chr_SSS.pvalues.ASM_lvl

# Remove the group effect
FL = np.repeat(np.nan, stData.shape[0])
for i in rnd_eff.keys():
    dt_dx = stData.patient_id == i
    if any(dt_dx):
        FL[dt_dx] = stData[rsp][dt_dx] - rnd_eff.get(i).Group
Data = stData.copy()
stData[rsp] = FL

plt.figure(figsize=[10,5])

# Plot all individuals: ASM level vs Seizure duration with MELM ---------------
plt.subplot (1,2,1)

# Plot the datapoints color codding seizure tyeps and ASM state they appear in
sns.scatterplot(data = stData, x="ASM_lvl", y=rsp, 
                hue = Cater, palette = palt, hue_order = Hue_Patt)

# Generate regresion for all datapoints
X = np.linspace(X_Lim_asm[0], X_Lim_asm[1])
Y = int_eff
Y += X * asm_eff * (asm_pvl < 0.05)
sns.lineplot(x=X, y=Y, color = "red", label = "ASM lvl. - Fixed eff.")

# Mark the normal-dose period
plt.fill_between(range(X_Lim_asm[0],X_Lim_asm[1]+1), Y_Lim2[0], Y_Lim2[1], 
                 np.array(range(X_Lim_asm[0],X_Lim_asm[1]+1))>=-1, 
                 alpha=0.1, color='k')

plt.xlim(X_Lim_asm)
plt.ylim(Y_Lim2)
plt.ylabel(Y_Label)
plt.xlabel("normalized ASM plasma concentration")

# Include sample size and ASM level fixed effect statistics
plt.title("All Seizures (" +
          str(stData.shape[0])+ " Szs. | "+
          str(len(np.unique(stData.patient_id)))+ " Pts.)" +
          "\n MELM - ASM lvl. eff=" + str(round(asm_eff,Stt_Dgt)) + 
          ", p=" + str(round(asm_pvl,Stt_Dgt)))

# Boxplot of seizure duration in patients with UND seizures -------------------
plt.subplot (1,2,2)

# Isolate individuals with any UND seizure
rm_dx = np.zeros(stData.shape[0])
for i in np.unique(stData.patient_id):
    idx = stData.patient_id == i
    sz_st = stData.SzSt[idx]
    if not any(sz_st == "UND"):
        rm_dx[idx] += 1
data = stData.drop(np.where(rm_dx)[0]).reset_index(drop=True)

# Alternative study: MELM to study Types and 
data_2 = Data.drop(np.where(rm_dx)[0]).reset_index(drop=True)
data_2.Sz_Class[data_2.Sz_Class == "UND_Under"] = "0_UND_Under"
data_2.Sz_Class[data_2.Sz_Class == "CMN_Under"] = "1_CMN_Under"
data_2.Sz_Class[data_2.Sz_Class == "CMN_Normo"] = "2_CMN_Normo"
mdl_Durt_TypAsm = smf.mixedlm(rsp + " ~ Sz_Class", data_2, groups=data_2.patient_id).fit()

# Boxplot of distribution and data poitns segregatin by seizure type and ASM perion
# Same color code and in previous plot
sns.boxplot(data = data, x=Cater, y=rsp, order = Axs_Patt,
            hue = Cater, palette = palt, hue_order = Hue_Patt,
            fill=False, flierprops={"marker": ""})

sns.stripplot(data = data, x=Cater, y=rsp, order = Axs_Patt,
              hue = Cater, palette = palt, hue_order = Hue_Patt,
              alpha = 0.7,jitter = 0.3)

xl = np.linspace(-0.5,2.5,7)
plt.fill_between(xl, Y_Lim2[0], Y_Lim2[1], xl>=1.5, alpha=0.1, color='k')
plt.xlim([xl[0],xl[-1]])

plt.ylim(Y_Lim2)

plt.ylabel("Duration(log10)")
plt.xlabel("")
plt.xticks([])
plt.title("Seizure Duration")

# Compare the difference between CMN and UND seizure duration
SSS_RankSum_UN = stats.ranksums (data[rsp][data.Sz_Class == "UND_Under"],
                              data[rsp][data.Sz_Class == "CMN_Under"])
N_UN = len(data[rsp][data.Sz_Class == "UND_Under"]) + len(data[rsp][data.Sz_Class == "CMN_Under"])
R_stat_UN = SSS_RankSum_UN.statistic/np.sqrt(N_UN)

# Compare the difference between CMN and UND seizure duration
SSS_RankSum_CM = stats.ranksums (data[rsp][data.Sz_Class == "UND_Under"],
                              data[rsp][data.Sz_Class == "CMN_Normo"])
N_CM = len(data[rsp][data.Sz_Class == "UND_Under"]) + len(data[rsp][data.Sz_Class == "CMN_Normo"])
R_stat_CM = SSS_RankSum_CM.statistic/np.sqrt(N_CM)

# Statistics of sieuzre type comparaison
plt.xlabel("UND vs CMN_Under - Rank Sum R=" + str(round(R_stat_UN,Stt_Dgt)) + ", p="+ str(round(SSS_RankSum_UN.pvalue,Stt_Dgt)) + '\n' +
           "UND vs CMN_Normo - Rank Sum R=" + str(round(R_stat_CM,Stt_Dgt)) + ", p="+ str(round(SSS_RankSum_CM.pvalue,Stt_Dgt)))

# Include sample size
plt.title("Individuals with UND seizures (" +
          str(data.shape[0])+ " Szs. | "+
          str(len(np.unique(data.patient_id)))+ " Pts.)")

if save_fig:
    plt.savefig(path_save_fig + "/RLT_ssDS_Durt_SzSt.svg", format='svg', 
                dpi=3000, bbox_inches='tight')
    plt.savefig(path_save_fig + "/PDF/RLT_ssDS_Durt_SzSt.pdf", format='pdf', 
                dpi=3000, bbox_inches='tight')

#%% Results: Seizures without Under-Dose Specific States  

# Full duration of patietns within Seizure States analysis
rsp = "DSS_Full"

Y_Lim1 = [0,3]
Y_Lim2 = [0,3]
Y_Label= "Duration (log10)"

X_Lim_asm = [-10,1]
X_Lim_tod = [0,24]

palt = "mako"
Hue_Patt = ["CMN_Normo","CMN_Under"]
Axs_Patt = ["CMN_Under","CMN_Normo"]
Cater = "Sz_Class"

# MELM for for duation against ASM levels and TofD
prd = "ASM_lvl + sTOD+cTOD + sTOD_h+cTOD_h"

# Remove all UND sieuzres
rm_dx = sz_Clust.SzSt == "UND"
# Clear dataset from non-relevant data
rm_dx_2 = np.logical_or(np.isnan(sz_Clust[rsp]),np.isinf(sz_Clust[rsp]))
rm_dx = np.logical_or(rm_dx,rm_dx_2)

stData = sz_Clust.drop(np.where(rm_dx)[0]).reset_index(drop=True)

# Fit the mixed effect model with patietn identifier group-level effect
mdl_Chr_CMN = smf.mixedlm(rsp + " ~ " + prd , stData, groups=stData.patient_id).fit()
mdl_Chr_CMN.summary()

# Relevent effect sizes and p-values for figure generation
rnd_eff = mdl_Chr_CMN.random_effects
int_eff = mdl_Chr_CMN.params.Intercept

asm_eff = mdl_Chr_CMN.params.ASM_lvl
asm_pvl = mdl_Chr_CMN.pvalues.ASM_lvl

# Remove the group effect
FL = np.repeat(np.nan, stData.shape[0])
for i in rnd_eff.keys():
    dt_dx = stData.patient_id == i
    if any(dt_dx):
        FL[dt_dx] = stData[rsp][dt_dx] - rnd_eff.get(i).Group
stData[rsp] = FL

plt.figure(figsize=[10,5])

# Plot CMN seizure duration vs ASM levels with MELM fixed effect --------------
plt.subplot(1,2,1)

# All data poitns collor codign the ASM period
sns.scatterplot(data = stData, x="ASM_lvl", y=rsp, 
                hue = Cater, palette = palt, hue_order = Hue_Patt)
# Regresion form MELM
X = np.linspace(X_Lim_asm[0], X_Lim_asm[1])
Y = int_eff
Y += X * asm_eff * (asm_pvl < 0.05)
sns.lineplot(x=X, y=Y, color = "red", label = "ASM lvl. - Fixed eff.")

# MArk the Mornal dose period
plt.fill_between(range(X_Lim_asm[0],X_Lim_asm[1]+1), Y_Lim2[0], Y_Lim2[1], 
                 np.array(range(X_Lim_asm[0],X_Lim_asm[1]+1))>=-1, 
                 alpha=0.1, color='k')

plt.xlim(X_Lim_asm)
plt.ylim(Y_Lim2)
plt.ylabel(Y_Label)
plt.xlabel("normalized ASM plasma concentration")

# OInclude sample size and MELM ASM level fixed-effect
plt.title("CMN Seizures (" +
          str(stData.shape[0])+ " Szs. | "+
          str(len(np.unique(stData.patient_id)))+ " Pts.)" +
          "\n MELM - ASM lvl. eff=" + str(round(asm_eff,Stt_Dgt)) + 
          ", p=" + str(round(asm_pvl,Stt_Dgt)))

# Boxplot of the difference between CMN seizures during Normal and Under dose -
plt.subplot(1,2,2)

# Boxplot of distribution and data poitns segregatin by seizure type and ASM perion
# Same color code and in previous plot
sns.stripplot(data = stData, x=Cater, y=rsp, order = Axs_Patt,
              hue = Cater, palette = palt, hue_order = Hue_Patt,
              alpha = 0.7, jitter = 0.3)

sns.boxplot(data = stData, x=Cater, y=rsp, order = Axs_Patt,
            hue = Cater, palette = palt, hue_order = Hue_Patt, 
            fill=False, flierprops={"marker": ""})

xl = np.array([-0.5,0,0.5,1,1.5])

# Mark dats corresponding to Normal-dose ASM period
plt.fill_between(xl, Y_Lim2[0], Y_Lim2[1], xl>=0.5, alpha=0.1, color='k')
plt.xlim([xl[0],xl[-1]])

plt.xlabel("")
plt.xticks([])

plt.ylim(Y_Lim2)
plt.ylabel(Y_Label)

# Statistical analysis of Normal-dose vs Under-dose CMN seizure duration
CMN_RankSum = stats.ranksums (stData[rsp][stData.ASM_stt == "Under"],
                              stData[rsp][stData.ASM_stt == "Normo"])
N = len(stData[rsp][stData.ASM_stt == "Under"]) + len(stData[rsp][stData.ASM_stt == "Normo"])
R_stat = CMN_RankSum.statistic/np.sqrt(N)

# Include statistical results
plt.xlabel("Under vs Normal - Rank Sum R=" + str(round(R_stat,Stt_Dgt)) + 
           ", p="+ str(round(CMN_RankSum.pvalue,Stt_Dgt)))

# Include sample size
plt.title("CMN Seizures (" +
          str(stData.shape[0])+ " Szs. | "+
          str(len(np.unique(stData.patient_id)))+ " Pts.)")

if save_fig:
    plt.savefig(path_save_fig + "/RLT_cmDS_Durt.svg", format='svg', 
                dpi=3000, bbox_inches='tight')
    plt.savefig(path_save_fig + "/PDF/RLT_cmDS_Durt.pdf", format='pdf', 
                dpi=3000, bbox_inches='tight')

#%% Supplementary: Under-Dose States Characterization
plt.figure(figsize=[10,5])
# Identify UND seizures
idx = sz_Clust.SzSt == "UND"

# Part 1: Sequence and relative duration of Under-dose Staes ------------------
plt.subplot(1,2,1)
ax = plt.gca()
ax.grid(axis='x',color='k', linestyle=':', linewidth=0.5)
xl = np.array([-20,0,20])
yl = [-0.5,2.5]

# Obtian distribution of Underdose/Common Sates on seizures
# 0.8 corresponds to non present state (n.p.)
pop_hist_seq = pd.DataFrame({"Seq":             ["n.p.", "1st", "2nd"],
                             "Under-dose Stt.": [-sum(sz_Clust.SqSS_Uniq[idx] == 0.8),
                                                 -sum(sz_Clust.SqSS_Uniq[idx] == 1),
                                                 -sum(sz_Clust.SqSS_Uniq[idx] == 2)],
                             "Common Stt.": [sum(sz_Clust.SqSS_Comn[idx] == 0.8),
                                             sum(sz_Clust.SqSS_Comn[idx] == 1),
                                             sum(sz_Clust.SqSS_Comn[idx] == 2)]})

# Barplots for Sequence distribution on each  ASM periods
sns.barplot(x='Under-dose Stt.', y='Seq', data=pop_hist_seq, order=["2nd", "1st", "n.p."], 
            orient='h', color="#EB6F3E", lw=0, label = 'Under-dose Stt.')

sns.barplot(x='Common Stt.', y='Seq', data=pop_hist_seq, order=["2nd", "1st", "n.p."], 
            orient='h', color="#531E77", lw=0, label = 'Common Stt.')

plt.legend()
plt.ylabel("sequence(index)")
plt.xlabel("counts(szs.)")

# Mark the Normal-dose related data
plt.fill_between(xl, yl[0], yl[1], xl>=0, alpha=0.1, color='k')

plt.xlim([xl[0],xl[-1]])
plt.ylim([yl[-1],yl[0]])

x_ticks = np.linspace(xl[0], xl[-1],int((xl[-1]-xl[0])/5+1))
plt.xticks(x_ticks,np.abs(x_ticks.astype(int)).astype(str))

plt.title("Seizure Stt. First Aperance")

# Part 2: Relative duration of Under-dose Staes -------------------------------
plt.subplot(1,2,2)

# Obtain the relative duration of Under-dose states
Data = 10**sz_Clust.DSS_Uniq / 10**sz_Clust.DSS_Full * 100

# Boxplot of distribution and data p[oints]
sns.boxplot(y=Data[idx],color = "#EB6F3E",fill=False, flierprops={"marker": ""})
sns.stripplot(y=Data[idx],alpha = 0.7,jitter = 0.3,color = "#EB6F3E", edgecolor='k')

plt.ylim([0,100])
plt.ylabel("Duration (%)")

# Median of the proporstion
plt.xlabel("median proportion = " +str(round(Data[idx].median(),3))+ "%")
plt.title("Under-dose states relative duration" )

if save_fig:
    plt.savefig(path_save_fig + "/SUP_UND_Seq_DurCont.svg", format='svg', 
                dpi=3000, bbox_inches='tight')
    plt.savefig(path_save_fig + "/PDF/SUP_UND_Seq_DurCont.pdf", format='pdf', 
                dpi=3000, bbox_inches='tight')

#%% Supplementary: Common States Characterization
plt.figure(figsize=[12,5])

# Part 1: Number of states within of Common States ----------------------------
plt.subplot(1,2,1)
ax = plt.gca()
ax.grid(axis='x',color='k', linestyle=':', linewidth=0.5)

# Identify CMN seizures with relevant data
rm_dx = sz_Clust.SzSt == "UND"
rm_dx_2 = np.logical_or(np.isnan(sz_Clust.DSS_Full),np.isinf(sz_Clust.DSS_Full))
rm_dx = np.logical_or(rm_dx,rm_dx_2)
Data = sz_Clust.drop(np.where(rm_dx)[0]).reset_index(drop = True)

# Obatin the distribution of the number of states of each seizure
nss = np.unique(Data.NSS_Comn)[::-1]
n_nss = np.zeros(nss.shape)
u_nss = np.zeros(nss.shape)
for i in range(0,len(nss)):
    n_nss[int(i)] = sum(Data.NSS_Comn[Data.ASM_stt == "Normo"] == nss[i])
    u_nss[int(i)] = sum(Data.NSS_Comn[Data.ASM_stt == "Under"] == nss[i])
pop_hist_nss = pd.DataFrame({"N. of States": nss,
                             "Under-dose N. Stt.":-u_nss,
                             "Normal-dose N. Stt.":n_nss})

# Barplots for Number of States distribution on each  ASM periods
sns.barplot(x='Under-dose N. Stt.', y='N. of States', data=pop_hist_nss,
            orient='h', color="#1B86AA", lw=0, label = 'Under-dose N. Stt.')

sns.barplot(x='Normal-dose N. Stt.', y='N. of States', data=pop_hist_nss,
            orient='h', color="#483B73", lw=0, label = 'Normal-dose N. Stt.')

plt.legend()
plt.ylabel("N. of States")
plt.xlabel("counts(szs.)")

xl = np.array([-45,0,45])
yl = [-0.5,len(nss)-0.5]
plt.fill_between(xl, yl[0], yl[1], xl>=0, alpha=0.1, color='k')

plt.xlim([xl[0],xl[-1]])
plt.ylim([yl[0],yl[-1]])

x_ticks = np.linspace(-40, 40,9)
plt.xticks(x_ticks,np.abs(x_ticks.astype(int)).astype(str))

plt.title("N. of States in CMN seizures(" +
          str(Data.shape[0])+ " Szs. | "+
          str(len(np.unique(Data.patient_id)))+ " Pts.)")

# Part 2: Contribution to length of states composing Common States ------------
plt.subplot(1,2,2)

# Identifu the contibution of each state to full duration
Full_Seg_L_CrCf = []
Full_Seg_L_CrpV = []

I = np.unique(Data.patient_id)
for Show_id in I:
    # Isolate patient data
    idx = np.where(Data.patient_id == Show_id)[0]
    SSS = [np.array(Data.sss_dx[i].split()).astype(int) for i in idx]
    
    # All states
    st_dx = np.unique(np.hstack(SSS))
    
    # Full duration of each seizure
    f_Lgt = np.array([len(s) for s in SSS])
    
    # Repository for Speramn Correlation
    cr_cf = np.zeros(len(st_dx)) #coefficient
    cr_pv = np.zeros(len(st_dx)) #pvalues
    cnt = 0
    for dx in st_dx:
        
        # Duration of the state on eahc seizure
        p_Lgt = np.array([sum(s==dx) for s in SSS])
        
        # Relationship between Full and Patial diration
        spearman = stats.spearmanr(f_Lgt,p_Lgt)
        cr_cf[cnt] = spearman.statistic
        # p-value as significance level [0:3]
        cr_pv[cnt] = sum([spearman.pvalue<0.05, 
                          spearman.pvalue<0.01,
                          spearman.pvalue<0.001])
        cnt += 1

    Full_Seg_L_CrCf.append(cr_cf)
    Full_Seg_L_CrpV.append(cr_pv)

# Scatterplot of Coefficients for each patietns color coding the signifiance
X = np.hstack([np.repeat(I[i], len(Full_Seg_L_CrCf[i])) for i in range(0,len(I))])
Y = np.hstack(Full_Seg_L_CrCf)
H = np.hstack(Full_Seg_L_CrpV)
sns.scatterplot(x = X.astype(str), y = Y, hue = H, 
                hue_norm = (0,4), palette="viridis_r")

plt.axhline(0,0,1,color='k', linestyle=":")

plt.title("Full Duration vs Common Sates Duration")
plt.ylabel("Spearman's Correlation Coefficient")

plt.legend(title='p-value', loc='lower left', labels=['','n.s.', '<0.05','<0.01','<0.001'])

x_ticks = ["id"+str(i) for i in I]
plt.xticks(np.linspace(0,len(I)-1,len(I)),x_ticks,rotation=30)

if save_fig:
    plt.savefig(path_save_fig + "/SUP_CMN_nSST_Corr.svg", format='svg', dpi=3000, bbox_inches='tight')
    plt.savefig(path_save_fig + "/PDF/SUP_CMN_nSST_Corr.pdf", format='pdf', dpi=3000, bbox_inches='tight')

#%% Supplementary: Effect of ASM types

# Get the full dataset seizure duration
rsp = "Duration_10"
rm_dx = np.logical_or(np.isnan(sz_Clust[rsp]),np.isinf(sz_Clust[rsp]))
stData = sz_Clust.drop(np.where(rm_dx)[0]).reset_index(drop=True)

plt.figure(figsize=[10,5])

palt = "mako"
Hue_Patt = ["CMN_Normo","CMN_Under","UND_Under"]

# ASM types on the dataset
# Type of ASM related to the target
# GB -> Gabba receptor
# SC -> Sodium Channels
# SV -> SVA1 receptor
# ML -> Multitarget
# OT -> Other medication (low representation)
drgt = np.array(['GB', 'ML', 'OT', 'SC', 'SV'])

# Usage of each ASM type across the dataset -----------------------------------
plt.subplot(1,2,1)

# Number of patietns with ASM type tapered
sbj_n = np.zeros(len(drgt))
for i in range(0,len(drgt)):
    idx = stData[drgt[i]] == 1
    sbj_n[i] = len(np.unique(stData.patient_id[idx]))
sns.barplot(x=drgt,y=sbj_n)
plt.title("ASM type used on tapering")

plt.ylim([0,len(np.unique(stData.patient_id))])
plt.ylabel("N. of Individuals")
plt.yticks(range(0, 29, 2))
plt.grid(visible=True,axis="y")

# Seizure duration segregates by ASM type -------------------------------------
plt.subplot(1,2,2)

for i in range(0,len(drgt)):
    # Identify patients with ASM type tapred
    idx = stData[drgt[i]] == 1
    
    # Distribution(boxplot) and datapoints fo all seiuzres/patietns identified
    # Color code seizure type
    sns.stripplot(x=i, y=stData[rsp][idx],legend=False, hue=stData.Sz_Class[idx], 
                  palette=palt, hue_order=Hue_Patt)
    sns.boxplot(x=i, y=stData[rsp][idx], fill=False, flierprops={"marker": ""}, color='k')

plt.xticks(range(0,len(drgt)),drgt)
plt.ylim([0,3.5])
plt.ylabel("Duration (log10)")
plt.title("Seizure duration by tapered ASM type")

if save_fig:
    plt.savefig(path_save_fig + "/SUP_ASM_Type.svg", format='svg', dpi=3000, bbox_inches='tight')
    plt.savefig(path_save_fig + "/PDF/SUP_ASM_Type.pdf", format='pdf', dpi=3000, bbox_inches='tight')

#%% Supplementary: Seizure dutation withouth Cyrc / Ultra dian effect
# Repica of 2nd resutls with changed MELM formula
# Fixed effects only related to ASM levels (No ToD variables included)
prd ="ASM_lvl"

rsp = "DSS_Full"

Y_Lim1 = [0,3]
Y_Lim2 = [0,3]
Y_Label= "Duration (log10)"

X_Lim_asm = [-10,1]
X_Lim_tod = [0,24]

palt = "mako"
Hue_Patt = ["CMN_Normo","CMN_Under","UND_Under"]
Axs_Patt = ["UND_Under","CMN_Under","CMN_Normo"]
Cater = "Sz_Class"

# Data of relevance
rm_dx = np.logical_or(np.isnan(sz_Clust[rsp]),np.isinf(sz_Clust[rsp]))
stData = sz_Clust.drop(np.where(rm_dx)[0]).reset_index(drop=True)

# Fit the MELM
mdl_Chr_Ams_nCrn = smf.mixedlm(rsp + " ~ " + prd , stData, groups=stData.patient_id).fit()
mdl_Chr_Ams_nCrn.summary()

rnd_eff = mdl_Chr_Ams_nCrn.random_effects
int_eff = mdl_Chr_Ams_nCrn.params.Intercept

asm_eff = mdl_Chr_Ams_nCrn.params.ASM_lvl
asm_pvl = mdl_Chr_Ams_nCrn.pvalues.ASM_lvl

# Remove Patient-level group effect
FL = np.repeat(np.nan, stData.shape[0])
for i in rnd_eff.keys():
    dt_dx = stData.patient_id == i
    if any(dt_dx):
        FL[dt_dx] = stData[rsp][dt_dx] - rnd_eff.get(i).Group
Data = stData.copy()
stData[rsp] = FL

# Plt seizure duration vs ASM level with MELM ASM level fixed effect
plt.figure(figsize=[5,5])

sns.scatterplot(data = stData, x="ASM_lvl", y=rsp, 
                hue = Cater, palette = palt, hue_order = Hue_Patt)
X = np.linspace(X_Lim_asm[0], X_Lim_asm[1])
Y = int_eff
Y += X * asm_eff * (asm_pvl < 0.05)
sns.lineplot(x=X, y=Y, color = "red", label = "ASM lvl. - Fixed eff.")
plt.fill_between(range(X_Lim_asm[0],X_Lim_asm[1]+1), Y_Lim2[0], Y_Lim2[1], 
                 np.array(range(X_Lim_asm[0],X_Lim_asm[1]+1))>=-1, 
                 alpha=0.1, color='k')
plt.xlim(X_Lim_asm)
plt.ylim(Y_Lim2)
plt.ylabel(Y_Label)
plt.xlabel("normalized ASM plasma concentration")
plt.title("All Seizures (" +
          str(stData.shape[0])+ " Szs. | "+
          str(len(np.unique(stData.patient_id)))+ " Pts.)" +
          "\n MELM - ASM lvl. eff=" + str(round(asm_eff,3)) + ", p=" + str(round(asm_pvl,3)))

if save_fig:
    plt.savefig(path_save_fig + "/SUP_Durt_SzSt_nCrn.svg", format='svg', 
                dpi=3000, bbox_inches='tight')
    plt.savefig(path_save_fig + "/PDF/SUP_Durt_SzSt_nCrn.pdf", format='pdf', 
                dpi=3000, bbox_inches='tight')



