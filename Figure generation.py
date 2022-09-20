#!/usr/bin/env python
# coding: utf-8

# # **Imported Modules**

# In[3]:


import matplotlib
import os
import numpy as np
import scipy.stats as stat
import pandas as pd
import seaborn as sns

from pandas import DataFrame as df
from matplotlib import pyplot as plt
from matplotlib.axes import Axes as ax


# # ***Figure 1***

# ## **Figure 1D. CRISPR screen result**

# In[7]:


positive = pd.read_csv("riger_pos_combine.csv")  
positive["i"]=positive.index
plot = sns.relplot(data=positive, x='i', y='P', aspect=1.5, 
                   hue='#CHROM', palette = 'Set3', size='size', legend=None, edgecolor=None) 
chrom_df=positive.groupby('#CHROM')['i'].median()
plot.ax.set_xticks(chrom_df)
plot.ax.set_xticklabels(chrom_df.index, rotation=90)
plot.ax.set_ylabel('Weighted average of p-value')
plot.ax.set_ylim(0,5.5)
plot.ax.axhline(y=2, linewidth = 1,linestyle="--",color="grey")
plot.fig.suptitle('Positive selection')
plt.savefig('positive selection.pdf', dpi=300)


# In[ ]:


negative = pd.read_csv("riger_neg_combine.csv")  
negative["i"]=negative.index
plot = sns.relplot(data=negative, x='i', y='P', aspect=1.5, 
                   hue='#CHROM', palette = 'Set3', size='size',legend=None, edgecolor=None) 
chrom_df=negative.groupby('#CHROM')['i'].median()
plot.ax.set_xticks(chrom_df)
plot.ax.set_xticklabels(chrom_df.index, rotation=90)
plot.ax.set_ylabel('Weighted average of p-value')
plot.ax.set_ylim(0,6)
plot.ax.axhline(y=2, linewidth = 1,linestyle="--",color="grey")
plot.fig.suptitle('Negative selection')
plot.fig.suptitle('Negative selection')
plt.savefig('negative selection.pdf', dpi=300)


# # ***Figure 2***

# ## **Figure 2A. ULK1 correlation with immune signatures**

# In[ ]:


rest={}
for i in data2.columns:
    rest[i]=[]


# In[ ]:


base="D:\\data\\"
for root, ds, fs in os.walk(base):
        for f in fs:
            print(base+f)
            file_list.append(f)
            data1=pd.read_csv(base+f)
            data1["Ensembl_ID"]=data1["Ensembl_ID"].apply(lambda x:x.split(".")[0])
            data1=data1.T
            data1.columns=data1.values[0]
            data1=data1[1:]
            data2=pd.read_csv("D:\\jupyter_workspace\\68ImmuneSigs.csv")
            data2=data2.set_index("Unnamed: 0")
            data_all=data1.join(data2,how="inner")
            t1=data_all["ENSG00000177169"]
            t2=data_all[data2.columns]
            for i in data2.columns:
                rest[i].append(stats.spearmanr(t1,t2[i], axis=1)[0])


# In[ ]:


s1=pd.DataFrame(rest)
s1.index=file_list
s1.to_csv("./ULK1_immune_sig.csv")


# In[ ]:


sns.clustermap(s1, cmap='RdBu_r',vmax=0.6, vmin=-0.6, row_cluster=False, col_cluster=False)
plt.title("ENSG00000177169 ULK1")
plt.show()


# ## **Figure 2C. ULK1 survival**

# In[ ]:


ulk1_metastasis=pd.read_csv('./ulk1_metastasis.csv', header=0)
ulk1_metastasis_col=ulk1_metastasis.T.iloc[0]
new_ulk1_metastasis=pd.DataFrame(ulk1_metastasis.T.iloc[1:].values,columns=ulk1_metastasis_col)
new_ulk1_metastasis["patient_sample"]=[i[0:16] for i in ulk1_metastasis.T.iloc[1:].index]


# In[ ]:


skcm_survival=pd.read_csv('./TCGA-SKCM.survival.tsv', sep='\t')
merge_ulk1=pd.merge(skcm_survival, new_ulk1_metastasis, on='patient_sample', how='inner')


# In[ ]:


font1 = {'family':'Arial', 'size' : 18}

km = KaplanMeierFitter()
T = merge_ulk1['OS.time'] 
E = merge_ulk1['OS']

high = (merge_ulk1['ULK1'] == 'high')
low = (merge_ulk1['ULK1'] == 'low')
ax = plt.subplot(111)
km.fit(T[high], event_observed=E[high], label="High expression of ULK1 (n=179)")
km.plot(ax=ax,ci_show=False, color='tomato')
km.fit(T[low], event_observed=E[low], label="Low expression of ULK1 (n=179)")
km.plot(ax=ax,ci_show=False, color='royalblue')
plt.title("TCGA-SKCM metastasis ULK1", font1)
plt.legend(prop=font2)
plt.savefig('TCGA-SKCM metastasis ULK1.pdf', dpi=300)


# ## **Figure 2D. IL10RB-DT correlation with immune signatures**

# In[ ]:


merge_skcm_sig=pd.read_csv('./skcm_dt.csv')
plt.figure(figsize=(8, 8))
sns.kdeplot(merge_skcm_sig['ENSG00000223799'].rank(method='average'),
        merge_skcm_sig['TAMsurr_TcClassII_ratio'].rank(method='average'),
        cmap='Blues', shade=True)
sns.regplot(
        x=merge_skcm_sig['ENSG00000223799'].rank(),
        y=merge_skcm_sig['TAMsurr_TcClassII_ratio'].rank(),
        scatter_kws={'s': 25, 'linewidths': .5, 'edgecolors': 'black',  'color': 'darkblue'},
        line_kws={'color': 'black', 'ls': '--', 'linewidth':1.5})
plt.xlim(0,len(merge_skcm_sig.index))
plt.ylim(0,len(merge_skcm_sig.index))
plt.xlabel("Rank in patients (IL10RB-DT expression)", font2)
plt.ylabel("Rank in patients (TAM signature score)", font2)
plt.title('SKCM', font2)
plt.savefig('DT-skcm-tam signature.pdf', dpi=600)
plt.show()


# In[ ]:


merge_thca=pd.read_csv("./thca-dt.csv")
plt.figure(figsize=(8, 8))
sns.kdeplot(
        merge_thca['ENSG00000223799.1'].rank(method='average'),
        merge_thca['CD8_CD68_ratio'].rank(method='average'),
        cmap='Blues', shade=True)
sns.regplot(
        x=merge_thca['ENSG00000223799.1'].rank(),
        y=merge_thca['CD8_CD68_ratio'].rank(),
        scatter_kws={'s': 25, 'linewidths': .5, 'edgecolors': 'black', 'color': 'darkblue'},
        line_kws={'color': 'black', 'ls': '--', 'linewidth':1.5})
plt.xlim(0,len(merge_thca.index))
plt.ylim(0,len(merge_thca.index))
plt.xlabel("Rank in patients (IL10RB-DT expression)", font2)
plt.ylabel("Rank in patients (CD8/CD68 signature score)", font2)
plt.title('THCA', font2)
plt.savefig('DT-thca-cd8cd68.pdf', dpi=600)
plt.show()


# In[ ]:


merge_laml=pd.read_csv('./laml_dt.csv')
plt.figure(figsize=(8, 8))
sns.kdeplot(
        merge_laml['ENSG00000223799.1'].rank(method='average'),
        merge_laml['CD8_CD68_ratio'].rank(method='average'),
        cmap='Blues', shade=True)
sns.regplot(
        x=merge_laml['ENSG00000223799.1'].rank(),
        y=merge_laml['CD8_CD68_ratio'].rank(),
        scatter_kws={'s': 35, 'linewidths': .5, 'edgecolors': 'black', 'color': 'darkblue'},
        line_kws={'color': 'black', 'ls': '--', 'linewidth':1.5})
plt.xlim(0,len(merge_laml.index))
plt.ylim(0,len(merge_laml.index))
plt.xlabel("Rank in patients (IL10RB-DT expression)", font2)
plt.ylabel("Rank in patients (CD8/CD68 signature score)", font2)
plt.title('LAML', font2)
plt.savefig('DT-LAML-cd8cd68.pdf', dpi=600)
plt.show()


# In[ ]:


merge_lgg=pd.read_csv('./lgg_dt.csv')
plt.figure(figsize=(8, 8))
sns.kdeplot(
        merge_lgg['ENSG00000223799.1'].rank(method='average'),
        merge_lgg['CD8_CD68_ratio'].rank(method='average'),
        cmap='Blues', shade=True)
sns.regplot(
        x=merge_lgg['ENSG00000223799.1'].rank(),
        y=merge_lgg['CD8_CD68_ratio'].rank(),
        scatter_kws={'s': 25, 'linewidths': .5, 'edgecolors': 'black', 'color': 'darkblue'},
        line_kws={'color': 'black', 'ls': '--', 'linewidth':1.5})
plt.xlim(0,len(merge_lgg.index))
plt.ylim(0,len(merge_lgg.index))
plt.xlabel("Rank in patients (IL10RB-DT expression)", font2)
plt.ylabel("Rank in patients (CD8/CD68 signature score)", font2)
plt.title('LGG', font2)
plt.savefig('DT-lgg-cd8cd68.pdf', dpi=600)
plt.show()


# ## **Figure 2F. IL10RB-DT expression and anti-PD1 treatment response**

# In[ ]:


treatment=pd.read_csv('./dt-pd1-response.csv')
fig, ax = plt.subplots(figsize=(1.5, 4))
sns.barplot(x="disease", y="expression", data=treatment, ci=68, capsize=.2, palette='Set2')
sns.swarmplot(x="disease", y="expression", data=treatment)
plt.savefig('dt-pd1-treatment.pdf', dpi=300)


# # ***Figure 3***

# ## **Figure 3D. volcano plot for IL10RB-DT activation**

# In[ ]:


plt.rcParams['font.sans-serif']='Arial'
pd_df = pd.read_csv('MEL526_IL10RBDT_Activation.csv')
result = pd.DataFrame()
result['x'] = pd_df['YW1_vs_YW5_Log2FC']
result['y'] = pd_df['p value']
x_threshold=0
y_threshold=1
result['group'] = 'black'
result.loc[(result.x > x_threshold)&(result.y > y_threshold),'group'] = 'tab:red' 
result.loc[(result.x < -x_threshold)&(result.y > y_threshold),'group'] = 'tab:blue'
result.loc[result.y < y_threshold,'group'] = 'grey' 

xmin=-3
xmax=3
ymin=0
ymax=4.2

fig = plt.figure(figsize=plt.figaspect(6/4)) 
ax = fig.add_subplot()
ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')
ax.scatter(result['x'], result['y'], s=2, c=result['group'])

plt.scatter([-0.6],[1.13],s=8,c='black')
plt.scatter([-0.76],[1.8],s=8,c='black')
plt.scatter([-0.55],[1.53],s=8,c='black')
plt.scatter([-0.63],[1.96],s=8,c='black')
plt.scatter([-0.71],[1.94],s=8,c='black')
plt.scatter([-0.69],[1.38],s=8,c='black')
plt.scatter([-0.53],[1.42],s=8,c='black')
plt.scatter([-0.72],[1.92],s=8,c='black')

plt.text(-0.6, 1.13, 'CXCL10', fontsize=13,  c='black', ha='right', va='bottom')
plt.text(-0.76, 1.7, 'IFITM1', fontsize=13,  c='black', ha='right', va='bottom')
plt.text(-0.55, 1.53, 'CXCL1', fontsize=13,  c='black', ha='right', va='bottom')
plt.text(-0.6, 1.96, 'IFI6', fontsize=13,  c='black', ha='left', va='bottom')
plt.text(-0.71, 2, 'CXCL3', fontsize=13, c='black', ha='right', va='bottom')
plt.text(-0.69, 1.38, 'ATF3', fontsize=13,  c='black', ha='right', va='bottom')
plt.text(-0.4, 1.42, 'PPP1R15A', fontsize=13, c='black', ha='left', va='bottom')
plt.text(-0.72, 1.85, 'CXCL2', fontsize=13, c='black', ha='right', va='bottom')

ax.set_ylabel('-log10(p value)')
ax.set_xlabel('log2 Fold Change')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.vlines(-x_threshold, ymin, ymax, color='dimgrey',linestyle='dashed', linewidth=1) 
ax.vlines(x_threshold, ymin, ymax, color='dimgrey',linestyle='dashed', linewidth=1) 
ax.hlines(y_threshold, xmin, xmax, color='dimgrey',linestyle='dashed', linewidth=1)

ax.set_xticks(range(-2, 3, 1)) 
ax.set_yticks(range(0,5,1))

plt.savefig('MEL526_DT_activation_0907.pdf', dpi=600)
plt.show()


# ## **Figure 3F. pathway analysis**

# In[ ]:


dt_table = pd.read_excel('dt_rnaseq.xlsx', sheet_name='NES')
dt_table['C'] = dt_table.groupby(['Term']).cumcount()+1
dt_table = dt_table.assign(new=dt_table.groupby(['Term'])['Term'].transform('count').astype('int'))
term = dt_table.apply(lambda row: row.Term if row.new==1 else row.Term + str(row.C), axis=1)
dt_table.drop(['C', 'new', 'Term'], inplace=True, axis=1)
dt_table.columns = [ 'MEL526 IL10RBDT-1', 'MEL526 IL10RBDT-2','MEL526 siIL10RBDT-1']

pvalue = pd.read_excel('dt_rnaseq.xlsx', sheet_name='FDR')
pvalue.drop(['Term'], inplace=True, axis=1)
pvalue.columns =  ['MEL526 IL10RBDT-1', 'MEL526 IL10RBDT-2','MEL526 siIL10RBDT-1']

dt_table.set_index([term], inplace=True)


# In[ ]:


xlabels = ['MEL526 IL10RBDT-1', 'MEL526 IL10RBDT-2','MEL526 siIL10RBDT-1']
ylabels = list(term)

fig, ax = plt.subplots(figsize=(4, 4), dpi=300,
                       facecolor='w', edgecolor='none')
divnorm = colors.TwoSlopeNorm(vmin=-2., vcenter=0., vmax=2)
im = ax.imshow(dt_table, cmap='coolwarm', norm=divnorm)

ax.set_xticks(np.arange(len(xlabels)))
ax.set_xticklabels(xlabels)
ax.set_yticks(np.arange(len(ylabels)))
ax.set_yticklabels(ylabels)

pmin = np.min(pvalue.to_numpy())
pmax = np.max(pvalue.to_numpy())
for i in range(len(ylabels)):
    for j in range(len(xlabels)):
        ax.scatter(j, i, s=200-((np.log10(
            pvalue.iloc[i, j]) - np.log10(pmin)) / (np.log10(pmax) - np.log10(pmin)))*180, c='black')
        ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, fill=False,
                           edgecolor='white', lw=1))

plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
         rotation_mode="anchor")

cax = fig.add_axes([ax.get_position().x1+0.185,
                   ax.get_position().y0, 0.02, ax.get_position().height])
clb = plt.colorbar(im, cax=cax)

clb.ax.set_title('NES', fontsize=8)

l1 = plt.scatter([], [], s=200-((np.log10(0.0001) - np.log10(pmin)) / (np.log10(pmax) - np.log10(pmin)))*180, color='black')
l2 = plt.scatter([], [], s=200-((np.log10(0.001) - np.log10(pmin)) / (np.log10(pmax) - np.log10(pmin)))*180, color='black')
l3 = plt.scatter([], [], s=200-((np.log10(0.01) - np.log10(pmin)) / (np.log10(pmax) - np.log10(pmin)))*180, color='black')
l4 = plt.scatter([], [], s=200-((np.log10(0.05) - np.log10(pmin)) / (np.log10(pmax) - np.log10(pmin)))*180, color='black')
labels = ['<0.0001', '<0.001', '<0.01', '<0.05']
fig.legend([l1, l2, l3, l4], labels, frameon=True, fontsize=12, title='FDR',
           borderaxespad=0, scatterpoints=1, bbox_to_anchor=(-0.5, -0.5))

fig.tight_layout()
fig.savefig('dt-RNAseq_pathway.pdf',dpi=300)
plt.show()


# # ***Figure 4***

# ## **Figure 4A. LINC01198 expression in normal and tumor tissues**

# In[ ]:


linc01198=pd.read_csv('./198-tcga-expression.csv')
fig, ax = plt.subplots(figsize=(5, 4))
sns.boxplot(x='type', y='ENSG00000231817.7', hue='tissue', data=linc01198, palette='Set2')
sns.stripplot(x="type", y="ENSG00000231817.7", hue="tissue", data=linc01198, palette="Set2", dodge=True)
plt.savefig('198-expression-tcga.pdf', dpi=300)


# # ***Figure 6***

# ## **Figure 6A. volcano plot for LINC01198 activation in MCF7**

# In[ ]:


plt.rcParams['font.sans-serif']='Arial'
pd_df = pd.read_csv('MCF7_LINC01198_Activation.csv')
result = pd.DataFrame()
result['x'] = pd_df['YW9_vs_YW10_Log2FC']
result['y'] = pd_df['p value']
x_threshold=0
y_threshold=1
result['group'] = 'black'
result.loc[(result.x > x_threshold)&(result.y > y_threshold),'group'] = 'tab:red' 
result.loc[(result.x < -x_threshold)&(result.y > y_threshold),'group'] = 'tab:blue' 
result.loc[result.y < y_threshold,'group'] = 'grey' 

xmin=-3
xmax=3
ymin=0
ymax=4.2

fig = plt.figure(figsize=plt.figaspect(6/4)) 
ax = fig.add_subplot()
ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')
ax.scatter(result['x'], result['y'], s=2, c=result['group'])

plt.scatter([1.65],[4],s=8,c='black')
plt.scatter([0.59],[1.46],s=8,c='black')
plt.scatter([0.69],[1.54],s=8,c='black')
plt.scatter([1.09],[3.07],s=8,c='black')
plt.scatter([0.75],[1.33],s=8,c='black')
plt.scatter([0.75],[1.33],s=8,c='black')
plt.scatter([0.81], [4], s=8, c='black')
plt.scatter([0.58], [2.22], s=8, c='black')
plt.scatter([0.54], [1.66], s=8, c='black')

ax.set_ylabel('-log10(p value)')
ax.set_xlabel('log2 Fold Change')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False) 

ax.vlines(-x_threshold, ymin, ymax, color='dimgrey',linestyle='dashed', linewidth=1)
ax.vlines(x_threshold, ymin, ymax, color='dimgrey',linestyle='dashed', linewidth=1) 
ax.hlines(y_threshold, xmin, xmax, color='dimgrey',linestyle='dashed', linewidth=1)

ax.set_xticks(range(-2, 3, 1)) 
ax.set_yticks(range(0,5,1)) 

plt.savefig('MCF7_LINC01198_activation_0907.pdf', dpi=600)
plt.show()


# ## **Figure 6B. pathway analysis for LINC01198 activation**

# In[ ]:


linc01198_table = pd.read_excel('198_rnaseq.xlsx', sheet_name='NES')
linc01198_table['C'] = linc01198_table.groupby(['Term']).cumcount()+1
linc01198_table = linc01198_table.assign(new=linc01198_table.groupby(['Term'])['Term'].transform('count').astype('int'))
term = linc01198_table.apply(lambda row: row.Term if row.new==1 else row.Term + str(row.C), axis=1)
linc01198_table.drop(['C', 'new', 'Term'], inplace=True, axis=1)
linc01198_table.columns = [ 'MEL526 LINC01198-1', 'MCF7 LINC01198-2','MEL526 siLINC01198-2']

pvalue = pd.read_excel('198_rnaseq.xlsx', sheet_name='FDR')
pvalue.drop(['Term'], inplace=True, axis=1)
pvalue.columns =  ['MEL526 LINC01198-1', 'MCF7 LINC01198-2','MEL526 siLINC01198-2']

linc01198_table.set_index([term], inplace=True)


# In[ ]:


xlabels = ['MEL526 LINC01198-1', 'MCF7 LINC01198-2','MEL526 siLINC01198-2']
ylabels = list(term)
fig, ax = plt.subplots(figsize=(4, 4), dpi=300,
                       facecolor='w', edgecolor='none')
divnorm = colors.TwoSlopeNorm(vmin=-2., vcenter=0., vmax=2)
im = ax.imshow(linc01198_table, cmap='coolwarm', norm=divnorm)

ax.set_xticks(np.arange(len(xlabels)))
ax.set_xticklabels(xlabels)
ax.set_yticks(np.arange(len(ylabels)))
ax.set_yticklabels(ylabels)

pmin = np.min(pvalue.to_numpy())
pmax = np.max(pvalue.to_numpy())
for i in range(len(ylabels)):
    for j in range(len(xlabels)):
        ax.scatter(j, i, s=200-((np.log10(
            pvalue.iloc[i, j]) - np.log10(pmin)) / (np.log10(pmax) - np.log10(pmin)))*180, c='black')
        ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, fill=False,
                           edgecolor='white', lw=1))

plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
         rotation_mode="anchor")

cax = fig.add_axes([ax.get_position().x1+0.185,
                   ax.get_position().y0, 0.02, ax.get_position().height])
clb = plt.colorbar(im, cax=cax)
clb.ax.set_title('NES', fontsize=8)
l1 = plt.scatter([], [], s=200-((np.log10(0.0001) - np.log10(pmin)) / (np.log10(pmax) - np.log10(pmin)))*180, color='black')
l2 = plt.scatter([], [], s=200-((np.log10(0.001) - np.log10(pmin)) / (np.log10(pmax) - np.log10(pmin)))*180, color='black')
l3 = plt.scatter([], [], s=200-((np.log10(0.01) - np.log10(pmin)) / (np.log10(pmax) - np.log10(pmin)))*180, color='black')
l4 = plt.scatter([], [], s=200-((np.log10(0.05) - np.log10(pmin)) / (np.log10(pmax) - np.log10(pmin)))*180, color='black')
labels = ['<0.0001', '<0.001', '<0.01', '<0.05']
fig.legend([l1, l2, l3, l4], labels, frameon=True, fontsize=12, title='FDR',
           borderaxespad=0, scatterpoints=1, bbox_to_anchor=(-0.5, -0.5))

fig.savefig('198-RNAseq0809_v2.pdf',dpi=300)
plt.show()


# # ***Figure S1***

# ## **Figure S1D. sgRNA frequencies for selected genes in different samples**

# In[ ]:


neg=pd.read_csv('./neg_selection_sgRNA.csv')
fig, ax = plt.subplots(figsize=(3, 5))
sns.boxplot(x='sample', y='log_count', data=neg, palette='Set2')
sns.stripplot(x="sample", y="log_count",  data=neg, palette="Set2", dodge=True)
plt.title('negative selection hits')
plt.savefig('negative selection hits sgRNA.pdf', dpi=300)


# In[ ]:


pos=pd.read_csv('./pos_selection_sgRNA.csv')
fig, ax = plt.subplots(figsize=(3, 5))
sns.boxplot(x='sample', y='log_count', data=pos, palette='Set2')
sns.stripplot(x="sample", y="log_count",  data=pos, palette="Set2", dodge=True)
plt.title('positive selection hits')
plt.savefig('positive selection hits sgRNA.pdf', dpi=300)


# # ***Figure S3***

# ## **Figure S3E. Immune-related genes regulated by IL10RB-DT**

# In[ ]:


df1 = pd.read_csv("dt_rnaseq_0323.csv",header=0, index_col=0)
plt.figure(figsize=(1,6))
ax= sns.heatmap(df1,cmap="bwr",center=0,vmin=-1.5,vmax=1.5,yticklabels=True)
ax.collections[0].colorbar.set_label("LFC",fontsize=10)

plt.savefig('dt-rnaseq-gene0323.pdf', dpi=300)

