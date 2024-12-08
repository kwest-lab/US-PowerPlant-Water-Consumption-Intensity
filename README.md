
# Water Consumption Intensity (WCI) Calculation for U.S. Power Plants

## Project Description  
This project calculates the **Water Consumption Intensity (WCI)** of electricity generation for U.S. power plants, using plant-level data provided by the Energy Information Administration (EIA). The code aggregates WCI values at the **Regional Transmission Organization (RTO)** level, allowing for insights into water use efficiency across regions.  

While this code is optimized for 2021 data, it can be adapted to analyze other years by replacing the input datasets with updated versions.

---

## Key Features  
- Computes water consumption intensity (gallons per megawatt-hour) for individual power plants.  
- Aggregates WCI values for major RTO regions across the U.S.  
- Handles various fuel types and energy sources, categorizing them into main groups (e.g., natural gas, petroleum).  
- Outputs intermediate datasets for further analysis and the final aggregated WCI values.  

---

## Requirements  
- **Python**: tested using Spyder script.  
- **Packages**:  
  - `pandas`  
  - `numpy`  

---

## Files  

### Input Data (Located in `input_data/`):  
1. **`data_CBG.csv`**  
   - Contains water consumption data for U.S. power plants (2021).  
   - Source: [EIA Electricity Data - Water](https://www.eia.gov/electricity/data/water/)  

2. **`EIA923_Schedules_2_3_4_5_M_12_2021_Final_Revision.csv`**  
   - Contains fuel consumption and electricity generation data for U.S. power plants (2021).  
   - Source: [EIA Form 923 Data](https://www.eia.gov/electricity/data/eia923/)  

3. **`fuel_codes_new.csv`**  
   - Categorizes fuel types into broader categories (e.g., natural gas, petroleum).  
   - Reference: [EIA Appendix C - Energy Source Grouping](https://www.eia.gov/electricity/monthly/pdf/AppendixC.pdf)  

4. **`RTO_PP_new.csv`**  
   - Assigns power plants to RTO regions based on GIS data processing.  
   - Sources:  
     - [Power Plant Locations - HIFLD](https://hifld-geoplatform.hub.arcgis.com/datasets/9dd630378fcf439999094a56c352670d_0/explore)  
     - [RTO Region Map - HIFLD](https://hifld-geoplatform.hub.arcgis.com/datasets/50f80920d36e435d9a34db2bd0fd3ad8_0/explore?location=35.610876%2C-95.679925%2C4.21)  

### Interim Data (Located in `interim_data/`):  
- Contains intermediate datasets generated during processing. Examples include WCI at the power plant level and energy production categorized by fuel type.

### Output Data (Located in `output_data/`):  
- Final aggregated WCI values by RTO regions.

---

## Usage  

1. Place all input data files in the `input_data/` directory.  
2. Open the Python script 'RTOs_new.py' in **Spyder** or your preferred IDE.  
3. Run the script to calculate WCI values and generate both interim and final outputs.  
4. Results are saved in the `output_data/` directory, with interim steps stored in `interim_data/`.  

---

## Citation  
If you use this code or the results generated from it, please cite the following sources for input data:  
- U.S. Energy Information Administration (EIA): [Electricity Data](https://www.eia.gov/)  
- Homeland Infrastructure Foundation-Level Data (HIFLD): [Geospatial Data](https://hifld-geoplatform.hub.arcgis.com/)  
