# JMONSEL-inverse

* **`AUTOMAP`:** the updated code of AUTOMAP for our problem.
* **`data`:** contains output of JMONSEL (`train_input`, `test_input`) and ground truth (`train_x_real`, `test_x_real`).
* **`fine_tuned_model/checkpoint`:** weights for fine-tuned AUTOMAP model with training part of JMONSEL output and ground truth images provided in `data` folder.
* **`pre_trained_model/checkpoint`:** weights for pre-trained AUTOMAP model with [MRI Training and Testing Data](https://www.dropbox.com/sh/fy5gnn6t1c6qgl2/AAAqIBMIaAlr4ZKLby-9u4QSa?dl=1).

## To Run
First you need to change configuration in `AUTOMAP/configs/inference_64x64_ex.json`. Give `loadmodel_dir` the path to the `fine_tuned_model/checkpoint` and add path for the `data` directory to `data_dir`.

Run the following command to create inferrence of the AUTOMAP output
```
python AUTOMAP/automap_main_inference.py -c AUTOMAP/configs/inference_64x64_ex.json
```

## JMONSEL-CODE
* More information about JMONSEL along with permission for it's code is given [JMONSEL-NIST](https://www.nist.gov/programs-projects/scanning-electron-microscopy)
* The folder jmonsel-code contains the output files of JMONSEL along with the code used to define the materials Silicon and Glassy Carbon. The code-file image-constructor.py was used to extract the information out of the textfiles generated as output of JMONSEL and converted it to a numpy matrix used for further processing.  
