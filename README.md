# perfusion_studies
Codes to Process the HP 13C Perfusion Studies

* The main code to process the 13C datasets is: process_c13_lungperfusion_20190807

* The path to the main folder for each dataset is saved as a structure in here: data_list_perfusion_studies_20190807
  * This makes it easier in future to run the analysis over multiple datasets if needed.
* Description for basic parameters
  * lb: line broadening in Hz
  * fit_spectra
    * if set to 0 then quicky phases the spectra and plots the time course for each metaboliate by calculating the area under the peak for the corresponding metabolite. 
    * If set to 1 then the code calls the "perfused_lung_13c_fit_20190807", which fits a Lorentizan functions to each peak to find the time-course. This options is slower but is more accurate. 
  * The results are saved as in "peakfit.txt" text file located the the corresponding dataset folder with multiple columns, which the code can load once run again. 
  * If the code is run for the second time, the results of the fit are loaded from the peakfit.txt file rather than applying the fit again. This saves time.
  

Codes to Process the HP 31P Perfusion Studies

* The main code to process the 31P dataset is: process_p31_lungperfusion_20190807

* This code for now is very preliminary and only loads the 31p datasets to show. We can edit this to quantify the peaks later.


