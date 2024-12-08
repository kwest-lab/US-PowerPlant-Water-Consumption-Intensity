import pandas as pd
import numpy as np
import os

# Get the current working directory
script_dir = os.getcwd()

# Define input and interim data directories
input_data_dir = os.path.join(script_dir, "input_data")
interim_data_dir = os.path.join(script_dir, "interim_data")
output_data_dir = os.path.join(script_dir, "output_data")

# Create interim_data and output_data directory if they don't exist
os.makedirs(interim_data_dir, exist_ok=True)
os.makedirs(output_data_dir, exist_ok=True)

#data preprocessing
def preprocess(x):
    return(x.replace(' ','0').replace(',',''))

#import cooling boiler generator data
dfW = pd.read_csv(os.path.join(input_data_dir, "data_CBG.csv"), parse_dates=False, low_memory=False, dtype=str)
dfW.columns = dfW.columns.str.replace(' ','_')
dfW['ID_SUP'] = dfW['State'] +'_'+dfW['Utility_ID']+'_'+dfW['Plant_Code']
headingsW = dfW.columns


#import EIA923 data
dfE = pd.read_csv(os.path.join(input_data_dir, "EIA923_Schedules_2_3_4_5_M_12_2021_Final_Revision.csv"), low_memory=False, dtype=str)
dfE.columns = dfE.columns.str.replace(' ','_')
dfE['ID_SUP'] = dfE['Plant_State'] +'_'+dfE['Operator_Id']+'_'+dfE['Plant_Id']
headingsE = dfE.columns

#Removing entries with 0 water withdrawals 
dfW['Water_Consumption_Volume_(Million_Gallons)'] =dfW['Water_Consumption_Volume_(Million_Gallons)'].astype(str)
dfW['Water_Consumption_Volume_(Million_Gallons)'] = pd.to_numeric(dfW['Water_Consumption_Volume_(Million_Gallons)'].apply(preprocess))
dfW['Water_Consumption_Volume_(Million_Gallons)'] = dfW['Water_Consumption_Volume_(Million_Gallons)'].astype(float)
dfW = dfW[dfW['Water_Consumption_Volume_(Million_Gallons)']>0]
dfW = dfW.reset_index(drop = True)  

#Removing entries with 0 electricity generation
dfE['Net_Generation_(Megawatthours)'] = pd.to_numeric(dfE['Net_Generation_(Megawatthours)'].apply(preprocess))
dfE['Net_Generation_(Megawatthours)'] = dfE['Net_Generation_(Megawatthours)'].astype(float) 
dfE = dfE[dfE['Net_Generation_(Megawatthours)']>0]

#Removing entries with 0 fuel consumption      
dfE['Elec_Fuel_Consumption_MMBtu'] = pd.to_numeric(dfE['Elec_Fuel_Consumption_MMBtu'].apply(preprocess))
dfE['Elec_Fuel_Consumption_MMBtu'] = dfE['Elec_Fuel_Consumption_MMBtu'].astype(float)        
dfE = dfE[dfE['Elec_Fuel_Consumption_MMBtu']>0]

#Filtering out the CHP plants
dfE = dfE[dfE['Sector_Name'].isin(['Electric Utility', 'NAICS-22 Non-Cogen'])]
#dfE = dfE[dfE['Sector_Name']=='Electric Utility' ]
dfE = dfE.reset_index(drop = True)

#Check common ID_SUP for both data sets
ID_SUP_W = dfW['ID_SUP'].unique().tolist()
ID_SUP_E = dfE['ID_SUP'].unique().tolist()

common = []
onlyW = []
onlyE = []
k=0 
while k<len(ID_SUP_W):
    if ID_SUP_W[k] in ID_SUP_E:
        common.append(ID_SUP_W[k])
        k=k+1 
    else:
        onlyW.append(ID_SUP_W[k])
        k=k+1
k=0 
while k<len(ID_SUP_E):
    if ID_SUP_E[k] in ID_SUP_W:
        k=k+1 
    else:
        onlyE.append(ID_SUP_E[k])
        k=k+1
        
        
#filtering data entries of cooling boiler generator data of common plants
k=0 
lstW = []
while k<(len(dfW['ID_SUP'])):
    if dfW['ID_SUP'][k] in common:
        new = dfW.iloc[k]
        lstW.append(new)
        k=k+1 
    else:
        k=k+1 
        
commonW = pd.DataFrame(lstW,columns=headingsW)
commonW = commonW.reset_index(drop = True)
commonW['Ftype'] = commonW['Generator_Primary_Energy_Source_Code']
       
 
#filtering data entries of eiadata data of common plants
k=0 
lstE = []
while k<(len(dfE['ID_SUP'])):
    if dfE['ID_SUP'][k] in common:
        new = dfE.iloc[k]
        lstE.append(new)
        k=k+1 
    else:
        k=k+1
         
commonE = pd.DataFrame(lstE,columns=headingsE)
commonE = commonE.reset_index(drop = True)
commonE['Ftype'] = commonE['Reported_Fuel_Type_Code']
    
#categorizing into main fueltypes water data
fuelcode = pd.read_csv(os.path.join(input_data_dir, "fuel_codes_new.csv"), dtype=str)
fuelcode.columns = fuelcode.columns.str.replace(' ','_')
fuellist = fuelcode['Energy_Source_Code'].tolist()

k=0 
while k<(len(commonW['Generator_Primary_Energy_Source_Code'])):
    if commonW['Generator_Primary_Energy_Source_Code'][k] in fuellist:
        commonW['Ftype'][k]= fuelcode['Fuel_Type'][fuellist.index(commonW['Generator_Primary_Energy_Source_Code'][k])]
        k=k+1 
    else:
        print('errorW')
        k=k+1 

#categorizing into main fueltypes electricity data
k=0 
while k<(len(commonE['Reported_Fuel_Type_Code'])):
    if commonE['Reported_Fuel_Type_Code'][k] in fuellist:
        commonE['Ftype'][k]= fuelcode['Fuel_Type'][fuellist.index(commonE['Reported_Fuel_Type_Code'][k])]
        k=k+1 
    else:
        print('errorE')
        k=k+1 


# annual water withdrawal values including ftype in ID
commonW['Unique'] = commonW['ID_SUP'] +'_'+commonW['Generator_ID']+'_'+commonW['Boiler_ID']+'_'+commonW['Cooling_ID']+'_'+commonW['Ftype']
commonW['Water_Consumption_Volume_(Million_Gallons)'] = commonW['Water_Consumption_Volume_(Million_Gallons)'].astype(float)        
a=np.shape(commonW['Unique'].unique())
Annual = np.zeros(((a[0]),2),dtype='<U64').astype(str)

k=0
Annual[:,0] = commonW['Unique'].unique()

for each in Annual[:,0]:
    Annual[k,1] = commonW.loc[commonW['Unique'] == each,'Water_Consumption_Volume_(Million_Gallons)'].sum()
    k=k+1 

Annual = pd.DataFrame(Annual) 
ID_split = Annual[0].str.split("_",expand =True)
Annual[2] = ID_split[6]
Annual[3] = ID_split[0]+'_'+ID_split[1]+'_'+ID_split[2]+'_'+ID_split[5]


# net electricity generation and fuel consumption including ftype in ID
commonE['Unique'] = commonE['ID_SUP'] +'_'+commonE['Ftype']
commonE['Net_Generation_(Megawatthours)'] = commonE['Net_Generation_(Megawatthours)'].astype(float)        
commonE['Elec_Fuel_Consumption_MMBtu'] = commonE['Elec_Fuel_Consumption_MMBtu'].astype(float)        

a1=np.shape(commonE['Unique'].unique())
netgen = np.zeros(((a1[0]),3),dtype='<U64').astype(str)
j=0

netgen[:,0] = commonE['Unique'].unique()

for each in netgen[:,0]:
    netgen[j,1] = commonE.loc[commonE['Unique'] == each,'Net_Generation_(Megawatthours)'].sum()
    netgen[j,2] = commonE.loc[commonE['Unique'] == each,'Elec_Fuel_Consumption_MMBtu'].sum()
    j=j+1 


# preprocessing annual data to identify plants with different fuel types with same cooling ID    
Annual_mod = Annual.drop(0,axis=1)
Annual_mod = Annual_mod.drop_duplicates()
Annual_mod = Annual_mod.sort_values(by=[3])
Annual_mod = Annual_mod.reset_index(drop = True)


#filtering plants that require percentage allocation and plants that doesnt require percentage allocation
perc_allc = []
nocalc =[]

k=0
while k<(len(Annual_mod[1])-1):
    if Annual_mod[3][k]== Annual_mod[3][k+1] and Annual_mod[2][k]!= Annual_mod[2][k+1]:
        perc_allc.append(Annual_mod.iloc[k])
        perc_allc.append(Annual_mod.iloc[k+1])
        k=k+1
    else:
        nocalc.append(Annual_mod.iloc[k])
        k=k+1 
        
#removing plants in nocalc which already in percentage calculation
Perc_Allc = pd.DataFrame(perc_allc)
Perc_Allc[0] = Perc_Allc[3]+'_'+Perc_Allc[2]   
No_calc = pd.DataFrame(nocalc)
No_calc[0] = No_calc[3]+'_'+No_calc[2] 
No_calc = No_calc.reset_index(drop=True)      
lstpercID = Perc_Allc[0].tolist()
k=0 
while k<(len(No_calc[0])):
    if No_calc[0][k] in lstpercID:
        No_calc = No_calc.drop(k)
        k=k+1 
    else:
        k=k+1 
No_calc = No_calc.reset_index(drop=True)         
    
# percentage allocation calculation of water withdrawals
Perc_Allc = Perc_Allc.drop_duplicates()
Perc_Allc[1] = Perc_Allc[1].astype(float)
Perc_Allc = Perc_Allc[Perc_Allc[1]>0]     
Perc_Allchead = Perc_Allc[3].str.split("_",expand = True)
Perc_Allc[4] = Perc_Allchead[0]+'_'+Perc_Allchead[1]+'_'+Perc_Allchead[2]+'_'+Perc_Allc[2]
Perc_Allc[5] = 0 #later fuel quantity 

Perc_Allc = Perc_Allc.reset_index(drop= True)
NetGen = pd.DataFrame(netgen)
listnet = NetGen[0].tolist()
k=0 
while k <(len(Perc_Allc[1])):
    if Perc_Allc[4][k] in listnet:
        Perc_Allc[5][k]= NetGen[2][listnet.index(Perc_Allc[4][k])]
        k=k+1 
    else:
        print('error')
        k=k+1 

a=np.shape(Perc_Allc)  
Perc_Allc[5] = Perc_Allc[5].astype(float)     
Perc_Allc1 = np.zeros(((a[0]),(a[1])),dtype='<U64').astype(str)
j=0

Perc_Allc1[:,0] = Perc_Allc[3] # 1st column plant ID
Perc_Allc1[:,1] = Perc_Allc[2] #2nd column fuel type
Perc_Allc1[:,2] = Perc_Allc[5].astype(float) #3rd column fuel quantity

for each in Perc_Allc1[:,0]:
    Perc_Allc1[j,3] = Perc_Allc.loc[Perc_Allc[3] == each,5].sum() #4th column total fuel quantity of plant
    j=j+1 
    
Perc_Allc1 =pd.DataFrame(Perc_Allc1)
Perc_Allc1[2] = Perc_Allc1[2].astype(float)
Perc_Allc1[3] = Perc_Allc1[3].astype(float)
Perc_Allc1[4] = Perc_Allc1[2]/Perc_Allc1[3] #5th column percentege per fuel type

Perc_Allc1[5] = Perc_Allc[1].astype(float) #6th column total water allocation
Perc_Allc1[6] = Perc_Allc1[4]*Perc_Allc1[5] #7th column percntage allocation per fuel type
 

#calculation of water intensities
Perc_Allchead = Perc_Allchead.reset_index(drop= True)
Perc_Allc1[7] = Perc_Allchead[0]+'_'+Perc_Allchead[1]+'_'+Perc_Allchead[2]+'_'+Perc_Allc[2]
Perc_Allc1[8] = 0 #later 9th column electricity generation 
k=0 
while k <(len(Perc_Allc1[7])):
    if Perc_Allc1[7][k] in listnet:
        Perc_Allc1[8][k]= NetGen[1][listnet.index(Perc_Allc1[7][k])]
        k=k+1 
    else:
        print('error')
        k=k+1 
Perc_Allc1[8] = Perc_Allc1[8].astype(float)
Perc_Allc1 = Perc_Allc1[Perc_Allc1[8]>0]
Perc_Allc1[9] = Perc_Allc1[6]/Perc_Allc1[8] #10th column water intensity
Perc_Allc1 = Perc_Allc1.fillna(0)


#data processing for calculation of water intensities that doesnt require percentage allocation
No_calc[1] = No_calc[1].astype(float)
No_calc = No_calc[No_calc[1]>0] 
No_calchead = No_calc[3].str.split("_",expand = True)
No_calc[4] = No_calchead[0]+'_'+No_calchead[1]+'_'+No_calchead[2]+'_'+No_calc[2]
No_calc[5] = -999 #later fuel quantity 
No_calc = No_calc.reset_index(drop= True)
k=0 
while k <(len(No_calc[4])):
    if No_calc[4][k] in listnet:
        No_calc[5][k]= NetGen[1][listnet.index(No_calc[4][k])]
        k=k+1 
    else:
        k=k+1 


#filtering  power plants with withdrawals yet no electricity generation
No_calc1 = []
nogeneration =[]

k=0
while k<(len(No_calc[5])):
    if No_calc[5][k]== -999:
        nogeneration.append(No_calc.iloc[k])
        k=k+1
    else:
        No_calc1.append(No_calc.iloc[k])
        k=k+1 


#calculation of water intensities that doesnt require percentage allocation
No_calc1 = pd.DataFrame(No_calc1)
No_calc1 = No_calc1.drop_duplicates()
No_calc1 = No_calc1.reset_index(drop= True)

No_calc1[5] = No_calc1[5].astype(float)
No_calc1 = No_calc1[No_calc1[5]>0]
No_calc1[6] = No_calc1[1]/No_calc1[5]


#combining both WI calculations
columns1 = [7,9]
columns2 = [4,6]
Perc_Allc_final = Perc_Allc1[columns1]
No_calc_final = No_calc1[columns2]
Perc_Allc_final.columns = ['Plant', 'WI']
No_calc_final.columns = ['Plant', 'WI']

WI = pd.concat([Perc_Allc_final,No_calc_final],axis=0)
WI = WI.reset_index(drop= True)


#power plant wise data for GRID map
WIGRID = WI.copy()


#adding electricity generation data
WIGRID['plant electricity'] = 0

k=0 
for each in WIGRID['Plant']:
    if WIGRID['Plant'][k] in listnet:
        WIGRID['plant electricity'][k] = NetGen[1][listnet.index(WIGRID['Plant'][k])]
        k=k+1 
    else:
        k=k+1 
WIGRID['ID_SUP'] = 0
lstWIGRID = WIGRID['Plant'].str.split("_",expand =True)
WIGRID['ID_SUP'] =lstWIGRID[0]+'_'+lstWIGRID[1]+'_'+lstWIGRID[2] 
WIGRID.index = WIGRID['ID_SUP']
lstWIGRID_ID = WIGRID['ID_SUP']

#adding RTO categories
dfRTO = pd.read_csv(os.path.join(input_data_dir, "RTO_PP_new.csv"), parse_dates=False, low_memory=False, dtype=str)
lstRTO = dfRTO['ID'].tolist()
WIGRID1 =[]
WIGRID1 = pd.DataFrame(WIGRID1)
WIGRID1['RTO'] = 0
k=0
j=0
for each in dfRTO['ID']:
    if each in lstWIGRID_ID:
        rtoplant = WIGRID.loc[[each]]
        WIGRID1 = pd.concat([WIGRID1,rtoplant],ignore_index=True)
        i= 1
        while i <= len(rtoplant):
            WIGRID1['RTO'][j] = dfRTO['RTO'][k]
            j=j+1
            i=i+1
        k=k+1
    else:
        k=k+1

lstWIGRID1 = WIGRID1['Plant'].str.split("_",expand =True)
WIGRID1['RTOcode'] = WIGRID1['RTO']+'_'+lstWIGRID1[3]

#calculation of RTO's average WI according to fuel type
f=np.shape(WIGRID1['RTOcode'].unique())
WIGRIDfinal = np.zeros(((f[0]),4),dtype='<U64').astype(str)
WIGRIDfinal[:,0] = WIGRID1['RTOcode'].unique()
WIGRID1['WI'] = WIGRID1['WI'].astype(float)
WIGRID1['plant electricity'] = WIGRID1['plant electricity'].astype(float)
WIGRID1['WI_2'] = WIGRID1['WI']*WIGRID1['plant electricity']

j=0
for each in WIGRIDfinal[:,0]:
    sumW = WIGRID1.loc[WIGRID1['RTOcode'] == each,'WI'].sum()
    sumW_2 = WIGRID1.loc[WIGRID1['RTOcode'] == each,'WI_2'].sum()
    ElecW_2 = WIGRID1.loc[WIGRID1['RTOcode'] == each,'plant electricity'].sum()
    countW = WIGRID1.loc[WIGRID1['RTOcode'] == each,'WI'].count()
    WIGRIDfinal[j,1] = sumW/countW
    WIGRIDfinal[j,2] = sumW_2/ElecW_2
    j=j+1 

# Filter rows from ISO-NE
iso_ne_row_NG = WIGRIDfinal[WIGRIDfinal[:, 0] == "ISO-NE_Natural Gas and Other Gases"]
iso_ne_row_P = WIGRIDfinal[WIGRIDfinal[:, 0] == "ISO-NE_Petroleum Products"]

# Create duplicates and modify GRID column
nyiso_row_NG = iso_ne_row_NG.copy()
nyiso_row_NG[:, 0] = "NYISO_Natural Gas and Other Gases"

nyiso_row_P = iso_ne_row_P.copy()
nyiso_row_P[:, 0] = "NYISO_Petroleum Products"

# Append the new rows to the original array
WIGRIDfinal = np.vstack((WIGRIDfinal, nyiso_row_NG, nyiso_row_P))

#calculation of total electricity generation according to fuel type


dfE['Ftype'] = 0
k=0 
while k<(len(dfE['Reported_Fuel_Type_Code'])):
    if dfE['Reported_Fuel_Type_Code'][k] in fuellist:
        dfE['Ftype'][k]= fuelcode['Fuel_Type'][fuellist.index(dfE['Reported_Fuel_Type_Code'][k])]
        k=k+1 
    else:
        print('errorE')
        k=k+1
dfE['Unique'] = dfE['ID_SUP'] +'_'+dfE['Ftype']
dfE['Net_Generation_(Megawatthours)'] = dfE['Net_Generation_(Megawatthours)'].astype(float)        

a1=np.shape(dfE['Unique'].unique())
netgen1 = np.zeros(((a1[0]),3),dtype='<U64').astype(str)
j=0

netgen1[:,0] = dfE['Unique'].unique()

for each in netgen1[:,0]:
    netgen1[j,1] = dfE.loc[dfE['Unique'] == each,'Net_Generation_(Megawatthours)'].sum()
    j=j+1 

#calculation of GRID's total electricity generation according to fuel type
NetGen2 = pd.DataFrame(netgen1)
NetGen2[1] = NetGen2[1].astype(float)
NetGen2 = NetGen2[NetGen2[1]>0]
lstNetGen2 = NetGen2[0].str.split("_",expand = True)
NetGen2[3] = lstNetGen2[0]+'_'+lstNetGen2[1]+'_'+lstNetGen2[2]
NetGen2[4] = 0
k=0
while k<(len(NetGen2[3])):
    if NetGen2[3][k] in lstRTO:
        NetGen2[4][k] = dfRTO['RTO'][lstRTO.index(NetGen2[3][k])]
        k=k+1 
    else:
        k=k+1 
NetGen2[4] = NetGen2[4].astype(str)    
NetGen2[5] = NetGen2[4]+'_'+lstNetGen2[3]
NetGen2[6] = lstNetGen2[2]
NetGen2[1] = NetGen2[1].astype(float)
NetGen2 = NetGen2[NetGen2[4]!='0']
j=0
for each in WIGRIDfinal[:,0]:
    WIGRIDfinal[j,3] = NetGen2.loc[NetGen2[5] == each,1].sum()
    j=j+1 

WIGRIDfinal = pd.DataFrame(WIGRIDfinal)
WIGRIDfinal[1] = WIGRIDfinal[1].astype(float) #WI_NA
WIGRIDfinal[2] = WIGRIDfinal[2].astype(float) #WI_WA
WIGRIDfinal[3] = WIGRIDfinal[3].astype(float) #Electricity
WIGRIDfinal[4] = WIGRIDfinal[1]*WIGRIDfinal[3] #WI_NA*electricity
WIGRIDfinal[5] = WIGRIDfinal[2]*WIGRIDfinal[3] #WI_WA*electricity

WIGRIDfinal.columns = ['GRIDcode','WI_NA','WI_WA,','electricity','WI_NA*electricity','WI_WA*electricity']
WIGRIDfinal = WIGRIDfinal[WIGRIDfinal['electricity']>0]
 
lstWIGRIDfinal = WIGRIDfinal['GRIDcode'].str.split("_",expand = True)
WIGRIDfinal['GRID'] = lstWIGRIDfinal[0]
g=np.shape(WIGRIDfinal['GRID'].unique())
WIGRIDfinal1 = np.zeros(((g[0]),5),dtype='<U64').astype(str)
WIGRIDfinal1[:,0] = WIGRIDfinal['GRID'].unique()

j=0
for each in WIGRIDfinal1[:,0]:
    WIGRIDfinal1[j,1] = WIGRIDfinal.loc[WIGRIDfinal['GRID'] == each,'WI_NA*electricity'].sum()
    WIGRIDfinal1[j,2] = WIGRIDfinal.loc[WIGRIDfinal['GRID'] == each,'WI_WA*electricity'].sum()
    j=j+1 

# Adding total electricity from all the energy sources
NetGen2[1] = NetGen2[1].astype(float)


j=0
for each in WIGRIDfinal1[:,0]:
    WIGRIDfinal1[j,3] = NetGen2.loc[NetGen2[4] == each,1].sum()
    j=j+1 

WIGRIDfinal1 = pd.DataFrame(WIGRIDfinal1)
WIGRIDfinal1[1] = WIGRIDfinal1[1].astype(float)
WIGRIDfinal1[2] = WIGRIDfinal1[2].astype(float)
WIGRIDfinal1[3] = WIGRIDfinal1[3].astype(float)
WIGRIDfinal1[4] = WIGRIDfinal1[1]/WIGRIDfinal1[3]
WIGRIDfinal1[5] = WIGRIDfinal1[2]/WIGRIDfinal1[3]
WIGRIDfinal1[6] = WIGRIDfinal1[4]*1000000
WIGRIDfinal1[7] = WIGRIDfinal1[5]*1000000


WIGRIDfinal1.columns = ['Balancing Authority','WI_NA*electricity','WI_WA*electricity','electricity','WI_NA in Mgal/MWh','WI_WA in Mgal/MWh','WI_NA in gal/MWh','WI_WA in gal/MWh']
NetGen.columns = ['Power plant','Net generation (MWh/year)','Electricity fuel consumption (MMBtu/year)']
Annual.columns = ['ID ','Water consumption (MGal/year)','Fuel type','Cooling ID']
NetGen2.columns = ['ID SUPF','Electricity generation (MWh/year)','','ID SUP','RTO','RTO and F type','Plant ID']

# Save processed datasets to interim_data
commonW.to_csv(os.path.join(interim_data_dir, "Common list of Water data.csv"), index=False)
commonE.to_csv(os.path.join(interim_data_dir, "Common list of Electricity data.csv"), index=False)
Annual.to_csv(os.path.join(interim_data_dir, "Annual water consumption of power plants.csv"), index=False)
NetGen.to_csv(os.path.join(interim_data_dir, "Annual electricity generation and fuel consumption of power plants.csv"), index=False)
Perc_Allc_final.to_csv(os.path.join(interim_data_dir, "Percentage allocation of power plats with different fuel types in same cooling ID.csv"), index=False)
NetGen2.to_csv(os.path.join(interim_data_dir, "Detailed annual electricity generation.csv"), index=False)
WI.to_csv(os.path.join(interim_data_dir, "Water consumption intensity_ Power plant level.csv"), index=False)


# Save porcesssed dataset to output_data
pd.DataFrame(WIGRIDfinal1).to_csv(os.path.join(output_data_dir, "WIGRIDfinal1.csv"), index=False)
