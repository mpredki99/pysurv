<details><summary style="font-size: 1.7em; font-weight: bold;">pysurv</summary>

---
(C) 2024, Michal Predki

**pysurv** is a Python package for adjusting surveying control networks.
The package supports importing data from *CSV* files and performing
**ordinary**, **weighted**, and **robust** least squares adjustment.
It also allows the **free adjustment approach** to be combined with
ordinary, weighted, and robust methods.
Additionally, you can mix these methods when adjusting
observations and reference points.
After completing the calculations, a detailed HTML report with
the adjustment results can be generated.

---

### Requirements

- **numpy**
- **scipy**
- **pandas**
- **matplotlib**

</details>

---

<details><summary style="font-size: 1.5em; font-weight: bold;">pysurv.customizations</summary>

---

[[source](src/pysurv/customizations.py)] This module defines default configurations for various parameters used in the adjustment of surveying control networks.

---

The following dictionaries are provided:

1. **methods**: Specifies default methods for adjusting observations and reference points in free adjustment.
    - **observations**: Default method for adjusting observations *(str)*.
    - **obs_c**: Tuning constant for robust adjustment of observations *(float, tuple of two floats, or None)*.
    - **free**: Default method for free adjustment *(str or None)*.
    - **free_c**: Tuning constant for robust free adjustment *(float, tuple of two floats, or None)*.

    If a robust method is provided and the tuning constant is None,
    a default theoretical value of the method will be used to calculate reweighting factors.

2. **default_measurement_sigma**: Contains default measurement sigma values for various quantities.
    - **sSD**: Default sigma value of 3D spatial distances *(tuple of two floats)*.
    - **sHD**: Default sigma value of 2D horizontal distances *(tuple of two floats)*.
    - **sVD**: Default sigma value of 1D vertical distances *(tuple of two floats)*.
    - **sdx**: Default sigma value of x-component of GNSS vectors *(tuple of two floats)*.
    - **sdy**: Default sigma value of y-component of GNSS vectors *(tuple of two floats)*.
    - **sdz**: Default sigma value of z-component of GNSS vectors *(tuple of two floats)*.
    - **sA**: Default sigma value of azimuthal angles *(float)*.
    - **sHZ**: Default sigma value of horizontal directions *(float)*.
    - **sVZ**: Default sigma value of vertical zenith angles *(float)*.
    - **sVH**: Default sigma value of vertical horizontal angles *(float)*.
    - **sP**: Default sigma value of the position of control points *(float)*.

3. **iterate_params**: Defines default parameters for iteration process in adjustment computations.
    - **tolerance**: Convergence tolerance *(float)*.
    - **max_iter**: Maximum number of iterations *(int)*.

4. **report_params**: Defines parameters used to create report with adjustment results.
    - **approx_precision**: number of decimals of approx coordinates and linear measurements *(int)*.
    - **adjusted_precision**: number of decimals of adjusted coordinates and linear measurements *(int)*.
    - **angle_unit**: unit to represent angles *(str)*.
    - **angle_precision**: number of decimals of angles *(int)*.

</details>

---

<details><summary style="font-size: 1.5em; font-weight: bold;">pysurv.Project</summary>

---

**Project**(path: *str*, methods: *str* = ***'default'***, measurement_errors: *str* = ***'default'***,
                 swap_xy: *bool* = ***False***, angle_unit: **tr = ***'grad'***,
     iterate_params: *str* = ***'default'***, report_params: *str* = ***'default'***)
---

[[source](src/pysurv/project_menager.py#L27)] A class to represent a surveying project that manages control points, measurements, and adjustment computations.

-----------------------------
**Initialization arguments:**

- **path**: *(str)*:

    The directory path to the project.

- **methods**: *(str, optional)*:

     Dictionary from customizations with the method set to use for adjustment.
     Defaults to 'default'.

- **measurement_errors**: *(str, optional)*:

    Dictionary from customizations with the default measurement sigma.
    Defaults to 'default'.

- **swap_xy**: *(bool, optional)*:

    Whether to swap the x and y coordinates when importing control points coordinates.
    Defaults to False.

- **angle_unit**: *(str, optional)*:

    The unit of the angles (e.g., 'grad', 'degree') in the measurements' dataset.
    Defaults to 'grad'.

- **iterate_params**: *(str, optional)*:

    Dictionary from customizations with parameters for the iterative process during
    adjustment computations. Defaults to 'default'.

- **report_params**: *(str, optional)*:

    Dictionary from customizations with parameters used for creating the adjustment
    results report. Defaults to 'default'.

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Attributes</summary>

- **path**: *(str)*:

    The directory path where project files are stored.

- **name**: *(str)*:

    The name of the project, derived from the project path.

- **controls**: *(Controls)*:

    The control points dataset.

- **measurements**: *(Measurements)*:

    The measurements' dataset.

- **methods**: *(dict)*:

    The methods to be used for adjustment.

- **default_measurement_errors**: *(dict)*:

    Default measurement errors, with angle errors converted to radians.

- **iterate_params**: *(dict)*:

    Parameters for the iterative process during adjustment computations.

- **report_params**: *(dict)*:

    Parameters used for creating the adjustment results report.

- **observation_equations**: *(dict or None)*:

    The matrices representing system of observation equations. At the time of initialization, it takes the value None.

- **adjustment_results**: *(dict or None)*:

    The results of the adjustment process. At the time of initialization, it takes the value None.

- **report**: *(Report or None)*:

    The report object containing the results of the adjustment process. At the time of initialization, it takes the value None.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Methods</summary>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">create_observation_system</summary>

---

**create_observation_system**()

---

[[source](src/pysurv/project_menager.py#L108)] Creates the system of observation equations.

Observations equations are represented as matrices.
This method updates the **observation_equations** attribute.

------------
**Returns:**

- **dict**: The system of observation equations:
  - **X**: *(np.ndarray)*: The coefficient matrix.
  - **Y**: *(np.ndarray)*: The vector of observed values.
  - **W**: *(np.ndarray, optional)*: The observations weight matrix.
  - **R**: *(np.ndarray, optional)*: The matrix of network inner constraints.
  - **sX**: *(np.ndarray, optional)*: The coefficient matrix for control points increments.
  - **sW**: *(np.ndarray, optional)*: The weight matrix for the control points.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjust</summary>

---

**adjust**(report_format: *str or None* = ***None***)

---

[[source](src/pysurv/project_menager.py#L133)] Performs the adjustment process on the dataset.

If the observation equations have not been created, they are generated first.
Adjustment computations are performed according to the selected methods.
After the calculation is completed, it is possible to generate a report
with the results to a file saved in the project folder.
This method updates the **adjustment_results** and **report** attributes.

--------------
**Arguments:**

- **report_format**: *(str, optional)*:

    The file format in which to export the report.
    If it is None, report will not be exported to the file.
    Defaults to None.

------------
**Returns:**

- **dict**: The results of the adjustment process:
  - **n_iter**: *(int)*: The number of iterations performed.
  - **sigma_zero**: *(list)*: List of sigma zero values from each iteration.
  - **pt_sigma_zero**: *(list)*: List of control points increments sigma zero values from each iteration.
  - **approx_coordinates**: *(Controls)*: The approximate control points coordinates before adjustment.
  - **adjusted_coordinates**: *(Controls)*: The adjusted control points coordinates.
  - **constraints_list**: *(list, optional)*: List of inner constraints applied to the network.

</details>

</details>

</details>

---

<details><summary style="font-size: 1.5em; font-weight: bold;">pysurv.modules</summary>

**pysurv** modules implement the tools needed to import data and perform least-squares network adjustment.

---

<details style="padding-left: 3%;"><summary style="font-size: 1.4em; font-weight: bold;">pysurv.modules.basic</summary>

[[source](src/pysurv/modules/basic.py)] This module provides utility functions for basic surveying calculations.

It includes functions to calculate the azimuth angle and to convert angles
between radians, degrees, and gradians (gons).

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Functions</summary>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">basic.azimuth</summary>

---

**azimuth**(first_point: *dict*, second_point: *dict*) -> **float**

---

[[source](src/pysurv/modules/basic.py#L16)] Calculate the azimuth angle between two points in a 2D plane, measured from the positive x-axis.

--------------
**Arguments:**

- **first_point**: *(dict)*:

    A dictionary with *'x'* and *'y'* coordinates for the first point.

- **second_point**: *(dict)*:

    A dictionary with *'x'* and *'y'* coordinates for the second point.

------------
**Returns:**

- **float**: The azimuth angle in radians, ranging from 0 to 2π.

-----------
**Raises:**

- **ValueError**: If the first and second points overlap.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">basic.to_rad</summary>

---

**to_rad**(angle: *float*, unit: *str* = ***'grad'***) -> **float**

---

[[source](src/pysurv/modules/basic.py#L47)] Convert an angle to radians from either gradians (gons) or degrees.

--------------
**Arguments:**

- **angle**: *(float)*:

    The angle to be converted.

- **unit**: *(str, optional)*:

    The unit of the angle *('grad', 'gon', or 'deg')*. Defaults to 'grad'.

------------
**Returns:**

- **float**: The angle in radians.

-----------
**Raises:**

- **ValueError**: If the unit is not 'grad', 'gon', or 'deg'.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">basic.from_rad</summary>

---

**from_rad**(angle: *float*, unit: *str* = ***'grad'***) -> **float**

---

[[source](src/pysurv/modules/basic.py#L75)] Convert an angle from radians to either gradians (gons) or degrees.

--------------
**Arguments:**

- **angle**: *(float)*:

    The angle to be converted.

- **unit**: *(str, optional)*:

    The unit to convert the angle to *('grad', 'gon', or 'deg')*. Defaults to 'grad'.

------------
**Returns:**

- **float**: The angle in the specified unit.

-----------
**Raises:**

- **ValueError**: If the unit is not 'grad', 'gon', or 'deg'.

</details>

</details>

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.4em; font-weight: bold;">pysurv.modules.importer</summary>

[[source](src/pysurv/modules/importer.py)] This module provides functionalities for import controls and measurements datasets.

Imported datasets are used to create instance of Controls and Measurements classes.

---

<details style="padding-left: 3%;"><summary style="font-size: 1.3em; font-weight: bold;">pysurv.modules.importer.CSV</summary>

[[source](src/pysurv/modules/importer.py#L18)] Class used to import datasets from CSV files.

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Methods</summary>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">controls</summary>

---

**controls**(path: *str*, swap_xy: *bool* = ***False***, *args, **kwargs) -> **Controls**

---

[[source](src/pysurv/modules/importer.py#L29)] Imports a CSV file containing control points coordinates and sigma values.

--------------
**Arguments:**

- **path**: *(str)*:

    Path to the CSV file containing the controls' dataset.

- **swap_xy**: *(bool)*:

    Whether to swap the values of *x* and *y* coordinates. Defaults to False.

- **\*args, \*\*kwargs**:

    Additional positional and keyword arguments passed to the pandas DataFrame initializer.

------------
**Returns:**

- **Controls**: An instance of the Controls class containing the controls dataset.

-----------
**Raises:**

- **ValueError**: If the provided path does not point to a valid file.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">measurements</summary>

---

**measurements**(path: *str*, angle_unit: *str* = ***'grad'***, *args, **kwargs) -> **Measurements**

---

[[source](src/pysurv/modules/importer.py#L58)] Imports a CSV file containing measurements and sigma values.

--------------
**Arguments:**

- **path**: *(str)*:

    Path to the CSV file containing the measurements' dataset.

- **angle_unit**: *(str)*:

    Unit of angular measurements in the dataset. Defaults to 'grad'.

- **\*args, \*\*kwargs**:

    Additional positional and keyword arguments passed to the pandas DataFrame initializer.

------------
**Returns:**

- **Measurements**: An instance of the Measurements class containing the measurements dataset.

-----------
**Raises:**

- **ValueError**: If the provided path does not point to a valid file.

</details>

</details>

</details>

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.4em; font-weight: bold;">pysurv.modules.adjustment</summary>

Adjustment submodule implements functions to define and solve a system of observation equations.

System of equations is solved by reducing them to a system of normal equations and solving this system by the least square method.

---

<details style="padding-left: 3%;"><summary style="font-size: 1.3em; font-weight: bold;">pysurv.modules.adjustment.computations</summary>

[[source](src/pysurv/modules/adjustment/computations.py)] This module contains functions to perform least squares adjustments for surveying control networks.

The adjustment process involves iterating over observation equations, updating weights
based on robust methods, and solving task with free adjustment approach.
The results include the adjusted coordinates, covariance matrices and other information
about adjustment process.

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Functions</summary>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.computations.adjust</summary>

---

**adjust**(controls: *Controls*, measurements: *Measurements*, matrices: *dict*,
           methods: *dict* = ***customizations.methods['default']***,
           iterate_params: *dict* = ***customizations.iterate_params['default']***) -> **dict**

---

[[source](src/pysurv/modules/adjustment/computations.py#L25)] Perform least squares adjustment on a set of control points and measurements.

--------------
**Arguments:**

- **controls**: *(Controls)*:

    Coordinates of control points.

- **measurements**: *(Measurements)*:

    The measurement dataset.

- **matrices**: *(dict)*:

    A dictionary containing matrices X, Y, and optionally W, R, sW, sX.

- **methods**: *(dict)*:

    Dictionary with methods to use to set observation and control point weights,
    and also with tuning constants. Defaults as 'default' in customisations.

- **iterate_params**: *(dict)*:

    Parameters for iteration process in computations. Defaults as 'default' in customisations.

------------
**Returns:**

- **dict**: A dictionary containing the adjustment results:
  - **n_iter**: *(int)*: The number of iterations performed.
  - **n_measurements**: *(int)*: Number of measurements.
  - **n_fixed_coords**: *(int)*: Number of fixed reference coordinates.
  - **n_sigma_coords**: *(int)*: Number of movable reference coordinates.
  - **n_unknowns**: *(int)*: Number of unknowns.
  - **r_norm**: *(int)*: Normalized residuals.
  - **b_norm**: *(int)*: Normalized increments.
  - **sigma_zero**: *(list)*: List of sigma zero values from each iteration.
  - **pt_sigma_zero**: *(list)*: List of control points increments sigma zero values from each iteration.
  - **approx_coordinates**: *(Controls)*: The approximate control points coordinates before adjustment.
  - **adjusted_coordinates**: *(Controls)*: The adjusted control points coordinates.
  - **constraints_list**: *(list)*: List of inner constraints applied to the network.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.computations.iterate</summary>

---

**iterate**(matrices: *dict*) -> **dict**

---

[[source](src/pysurv/modules/adjustment/computations.py#L146)] Perform iteration for least squares adjustment.

Solves a system of equations and determines covariance matrices.

--------------
**Arguments:**

- **matrices**: *(dict)*:

    A dictionary containing matrices:

  - **X**: *(np.ndarray)*: The coefficient matrix.
  - **Y**: *(np.ndarray)*: The vector of observed values.
  - **W**: *(np.ndarray, optional)*: The observations weight matrix.
  - **R**: *(np.ndarray, optional)*: The matrix of network inner constraints.
  - **sX**: *(np.ndarray, optional)*: The coefficient matrix for control points increments.
  - **sW**: *(np.ndarray, optional)*: The weight matrix for the control points.

------------
**Returns:**

- **dict**: A dictionary containing the iteration results:
  - **increments**: *(np.ndarray)*: Vector of calculated increments to coordinates.
  - **obs_residuals**: *(np.ndarray)*: Vector of observation residuals.
  - **sigma_zero_sq**: *(np.ndarray)*: Residual variance.
  - **N**: *(np.ndarray)*: Matrix of normal equations.
  - **L**: *(np.ndarray)*: Vector of dependent variable.
  - **cov_b**: *(np.ndarray)*: Covariance matrix of increments.
  - **cov_Y**: *(np.ndarray)*: Covariance matrix of observations.
  - **cov_r**: *(np.ndarray)*: Covariance matrix of residuals.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.computations.normal_equations</summary>

---

**normal_equations**(X: *np.ndarray*, W: *np.ndarray*, Y: *np.ndarray*,
                     sW: *np.ndarray*, R: *np.ndarray*, sX: *np.ndarray*) -> **Tuple[np.ndarray, np.ndarray]**

---

[[source](src/pysurv/modules/adjustment/computations.py#L211)] Determine the normal equations and dependent variable for the adjustment process.

--------------
**Arguments:**

- **X**: *(np.ndarray)*:
  
    Coefficient matrix.

- **W**: *(np.ndarray)*:

    The observations weight matrix.

- **Y**: *(np.ndarray)*:

    The vector of observed values.

- **sW**: *(np.ndarray)*:

    The weight matrix for the control points.

- **R**: *(np.ndarray)*:

    The matrix of network inner constraints.

- **sX**: *(np.ndarray)*:

    The coefficient matrix for control points increments.

------------
**Returns:**

- Tuple containing normal equations and dependent variable:
  - **np.ndarray**: The matrix of normal equations.
  - **np.ndarray**: The vector of dependent variable.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.computations.update_weights</summary>

---

**update_weights**(matrices: *dict*, iterate_results: *dict*, methods: *dict*, pt_sigma_zero_sq: *float*) -> **None**

---

[[source](src/pysurv/modules/adjustment/computations.py#L248)] Update the weight matrices based on robust methods during the iteration process.

This function modifies the input array in-place without returning a value.

--------------
**Arguments:**

- **matrices**: *(dict)*:

    A dictionary containing weight matrices.

- **iterate_results**: *(dict)*:

    The results from the current iteration.

- **methods**: *(dict)*:

    Dictionary specifying the robust methods and tuning constants to use.

- **pt_sigma_zero_sq**: *(float)*:

    Control points' increments squared sigma zero parameter.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.computations.reweight</summary>

---

**reweight**(weights: *np.ndarray*, residuals: *np.ndarray*, residuals_variances: *np.ndarray*, method: *str*, kwargs: *dict*) -> **None**

---

[[source](src/pysurv/modules/adjustment/computations.py#L302)] Recalculate the weights based on the normalized residuals and a selected robust method.

This function modifies the input array in-place without returning a value.

--------------
**Arguments:**

- **weights**: *(np.ndarray)*:

    The current weight matrix.

- **residuals**: *(np.ndarray)*:

    The residual values to be reweighted.

- **residuals_variances**: *(np.ndarray)*:

    The variances of the residuals.

- **method**: *(str)*:

    The robust method to use for reweighting.

- **\*\*kwargs**: *(dict)*:

    Keyword arguments required by the robust method.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.computations.normalize_residuals</summary>

---

**normalize_residuals**(residuals: *np.ndarray*, residuals_variances: *np.ndarray*) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/computations.py#L322)] Normalize values of residuals based on their values and variances.

--------------
**Arguments:**

- **residuals**: *(np.ndarray)*:

    The residual values.

- **residuals_variances**: *(np.ndarray)*:

    The variances of the residuals.

------------
**Returns:**

- **np.ndarray**: Array of normalized values of residuals.

</details>

</details>

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.3em; font-weight: bold;">pysurv.modules.adjustment.matrices</summary>

[[source](src/pysurv/modules/adjustment/matrices.py)] A module for creating and managing matrices used in the adjustment of surveying control networks.

It supports the creation of observation equation matrices, weight matrices, and constraint matrices
for free adjustment, as well as matrices for adjustment including the standard deviations of control points.

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Functions</summary>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.matrices.equations_system</summary>

---

**equations_system**(controls: *Controls*, measurements: *Measurements*,
                     methods: *dict* = ***customizations.methods['default']***,
                     default_sigma: *dict* = ***customizations.default_measurement_sigma['default']***) -> **dict**

---

[[source](src/pysurv/modules/adjustment/matrices.py#L23)] Creates a system of observation equations for adjustment of surveying control network.

System of equations is represented as the matrices.

--------------
**Arguments:**

- **controls**: *(Controls)*:

    The control points dataset.

- **measurements**: *(Measurements)*:

    The measurement dataset.

- **methods**: *(dict, optional)*:

    Dictionary with methods to use to set observation and control point weights, and also with tuning constants.
    Defaults as 'default' in customisations.

- **default_sigma**: *(dict, optional)*:

    Default standard deviations for measurements and point locations. Defaults as 'default' in customisations.

------------
**Returns:**

- **dict**: A dictionary containing matrices:
  - **X**: *(np.ndarray)*: The coefficient matrix.
  - **Y**: *(np.ndarray)*: The vector of observed values.
  - **W**: *(np.ndarray, optional)*: The observations weight matrix.
  - **R**: *(np.ndarray, optional)*: The matrix of network inner constraints.
  - **sX**: *(np.ndarray, optional)*: The coefficient matrix for control points increments.
  - **sW**: *(np.ndarray, optional)*: The weight matrix for the control points.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.matrices.XYW_matrices</summary>

---

**XYW_matrices**(controls: *Controls*, measurements: *Measurements*, calculate_W: *bool* = ***True***,
                 default_sigma: *dict* = ***customizations.default_measurement_sigma['default']***) -> **dict**

---

[[source](src/pysurv/modules/adjustment/matrices.py#L86)] Creates a system of observation equations for surveying measurements.

This function iterates through each pair of points specified in the measurements' dataset,
calculates the necessary observation equations using the appropriate observation function,
and builds the X (coefficient matrix), Y (observation vector), and optionally W (diagonal weight matrix) matrices.

--------------
**Arguments:**

- **controls**: *(Controls)*:

    Coordinates of control points.

- **measurements**: *(Measurements)*:

    The measurement dataset.

- **calculate_W**: *(bool, optional)*:

    If True, the function calculates the weight matrix (W) for the observations. Defaults to True.

- **default_sigma**: *(dict, optional)*:

    Default standard deviations for measurements and point locations. Defaults as 'default' in customisations.

------------
**Returns:**

- **dict**: A dictionary containing matrices:
  - **X**: *(np.ndarray)*: The coefficient matrix.
  - **Y**: *(np.ndarray)*: The vector of observed values.
  - **W**: *(np.ndarray, optional)*: The observations weight matrix.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.matrices.apply_inner_constraints</summary>

---

**apply_inner_constraints**(controls: *Controls*, measurement_types: *pd.Index*) -> **dict**

---

[[source](src/pysurv/modules/adjustment/matrices.py#L162)] Applies inner constraints for free network adjustment.

Function creates the R matrix (matrix of inner constraints).

--------------
**Arguments:**

- **controls**: *(Controls)*:

    Coordinates of control points.

- **measurement_types**: *(pd.Index)*:

    A list of measurement types in dataset.

------------
**Returns:**

- **dict**: A dictionary with 'R' matrix and list of applied constraints:
  - **R**: *(np.ndarray)*: Matrix of inner constraints.
  - **constraints**: *(list)*: List of applied constraints.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.matrices.calculate_sW</summary>

---

**calculate_sW**(controls: *Controls*, calculate_sX: *bool* = ***True***,
                 default_pt_sigma: *float* = ***customizations.default_measurement_sigma['default']['sP']***) -> **dict**

---

[[source](src/pysurv/modules/adjustment/matrices.py#L239)] Determines the sW and optionally sX matrices.

Calculates the diagonal matrix of control points' weights (sW) and, optionally,
the coefficient matrix for increments (sX) for adjustments that include control points' standard deviations.

--------------
**Arguments:**

- **controls**: *(Controls)*:

    Coordinates of control points and their standard deviations.

- **calculate_sX**: *(bool)*:

    If True, the function calculates the sX matrix for coefficient matrix increments. Defaults to True.

- **default_pt_sigma**: *(float, optional)*:

    Default standard deviation of points location. Defaults as 'default sP' in customisations.

------------
**Returns:**

- **dict**: A dictionary containing matrices:
  - **sW** *(np.ndarray)*: The weight matrix for the control points.
  - **sX** *(np.ndarray, optional)*: The coefficient matrix for control points increments.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.matrices.approx_orientation</summary>

---

**approx_orientation**(controls: *Controls*, from_to_hz: *pd.Series*) -> *pd.Series*

---

[[source](src/pysurv/modules/adjustment/matrices.py#L295)] Computes the approximate orientation constant at control points.

--------------
**Arguments:**

- **controls**: *(Controls)*:

    Coordinates of control points.

- **from_to_hz**: *(pd.Series)*:

    A Series with MultiIndex of ('FROM', 'TO') and values of horizontal directions.

------------
**Returns:**

- **pd.Series**: A Series containing the computed orientation constant for each point.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.matrices.set_obs_weight</summary>

---

**set_obs_weight**(row: *pd.Series*, measurement_type: *str*, measurement_value: *float*,
                   default_sigma: *dict* = ***customizations.default_measurement_sigma['default']***) -> **float**

---

[[source](src/pysurv/modules/adjustment/matrices.py#L321)] Computes the observation weight based on the measurement type and its standard deviation.

If an estimated sigma value is given in the dataset, the point weight is calculated based on it.
If this value is not present, the point weight is calculated based on the default value defined in
the default_sigma dictionary.

--------------
**Arguments:**

- **row**: *(pd.Series)*:

    A row from a measurements' dataset.

- **measurement_type**: *(str)*:

    The type of measurement.

- **measurement_value**: *(float)*:

    The value of the measurement.

- **default_sigma**: *(dict)*:

    Default standard deviations for measurements and points' localization.

------------
**Returns:**

- **float**: The observation weight.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.matrices.inv_sq_error</summary>

---

**inv_sq_error**(sigma: *float*) -> **float**

---

[[source](src/pysurv/modules/adjustment/matrices.py#L352)] Sets the observation weight based on the standard deviation.

The weight is defined as the inverse of the square of the standard deviation.

--------------
**Arguments:**

- **sigma**: *(float)*:

    Standard deviation value.

------------
**Returns:**

- **float**: Weight.

</details>

</details>

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.3em; font-weight: bold;">pysurv.modules.adjustment.observation_equations</summary>

[[source](src/pysurv/modules/adjustment/observation_equations.py)] The module contains functions for determining the coefficients of the observation equations of surveying measurements.

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Functions</summary>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.observation_equations.SD_obs_eq</summary>

---

**SD_obs_eq**(meas_SD: *float*, idx_from: *dict*, idx_to: *dict*, coord_differences: *dict*, X_row: *np.ndarray*) -> **Tuple[np.ndarray, np.ndarray]**

---

[[source](src/pysurv/modules/adjustment/observation_equations.py#L16)] Calculates the observation equation for slope distance *(SD)*.

--------------
**Arguments:**

- **meas_SD**: *(float)*:

    Measured slope distance.

- **idx_from**: *(dict)*:

    Indices for the station point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the station.
  - **y**: *(int)*: The index of the y-coordinate increment column of the station.
  - **z**: *(int)*: The index of the z-coordinate increment column of the station.

- **idx_to**: *(dict)*:

    Indices for the aim point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the aim point.
  - **y**: *(int)*: The index of the y-coordinate increment column of the aim point.
  - **z**: *(int)*: The index of the z-coordinate increment column of the aim point.

- **coord_differences**: *(dict)*:

    Coordinate differences between the points:

  - **dx**: *(float)*: The difference of the x-coordinates of the station and the aim point.
  - **dy**: *(float)*: The difference of the y-coordinates of the station and the aim point.
  - **dz**: *(float)*: The difference of the z-coordinates of the station and the aim point.

- **X_row**: *(np.ndarray)*:

    The X matrix row to be populated.

------------
**Returns:**

- **tuple**: A tuple of arrays:
  - *(np.ndarray)*: updated X_row.
  - *(np.ndarray)*: updated Y_row.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.observation_equations.HD_obs_eq</summary>

---

**HD_obs_eq**(meas_HD: *float*, idx_from: *dict*, idx_to: *dict*, coord_differences: *dict*, X_row: *np.ndarray*) -> **Tuple[np.ndarray, np.ndarray]**

---

[[source](src/pysurv/modules/adjustment/observation_equations.py#L64)] Calculates the observation equation for horizontal distance *(HD)*.

--------------
**Arguments:**

- **meas_HD**: *(float)*:

    Measured horizontal distance.

- **idx_from**: *(dict)*:

    Indices for the station point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the station.
  - **y**: *(int)*: The index of the y-coordinate increment column of the station.

- **idx_to**: *(dict)*:

    Indices for the aim point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the aim point.
  - **y**: *(int)*: The index of the y-coordinate increment column of the aim point.

- **coord_differences**: *(dict)*:

    Coordinate differences between the points:

  - **dx**: *(float)*: The difference of the x-coordinates of the station and the aim point.
  - **dy**: *(float)*: The difference of the y-coordinates of the station and the aim point.

- **X_row**: *(np.ndarray)*:

    The X matrix row to be populated.

------------
**Returns:**

- **tuple**: A tuple of arrays:
  - *(np.ndarray)*: updated X_row.
  - *(np.ndarray)*: updated Y_row.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.observation_equations.dx_obs_eq</summary>

---

**dx_obs_eq**(meas_dx: *float*,
              idx_from: *dict*,
              idx_to: *dict*,
              coord_differences: *dict*,
              X_row: *np.ndarray*) -> **Tuple[np.ndarray, np.ndarray]**

---

[[source](src/pysurv/modules/adjustment/observation_equations.py#L107)] Calculates the observation equation for the x-coordinate component of the GNSS vector *(dx)*.

--------------
**Arguments:**

- **meas_dx**: *(float)*:

    Measured x-coordinate difference.

- **idx_from**: *(dict)*:

    Indices for the station point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the station.

- **idx_to**: *(dict)*:

    Indices for the aim point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the aim point.

- **coord_differences**: *(dict)*:

    Coordinate differences between the points:

  - **dx**: *(float)*: The difference of the x-coordinates of the station and the aim point.

- **X_row**: *(np.ndarray)*:

    The X matrix row to be populated.

------------
**Returns:**

- **tuple**: A tuple of arrays:
  - *(np.ndarray)*: updated X_row.
  - *(np.ndarray)*: updated Y_row.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.observation_equations.dy_obs_eq</summary>

---

**dy_obs_eq**(meas_dy: *float*,
              idx_from: *dict*,
              idx_to: *dict*,
              coord_differences: *dict*,
              X_row: *np.ndarray*) -> **Tuple[np.ndarray, np.ndarray]**

---

[[source](src/pysurv/modules/adjustment/observation_equations.py#L143)] Calculates the observation equation for the y-coordinate component of the GNSS vector *(dy)*.

--------------
**Arguments:**

- **meas_dy**: *(float)*:

    Measured y-coordinate difference.

- **idx_from**: *(dict)*:

    Indices for the station point:

  - **y**: *(int)*: The index of the y-coordinate increment column of the station.

- **idx_to**: *(dict)*:

    Indices for the aim point:

  - **y**: *(int)*: The index of the y-coordinate increment column of the aim point.

- **coord_differences**: *(dict)*:

    Coordinate differences between the points:

  - **dy**: *(float)*: The difference of the y-coordinates of the station and the aim point.

- **X_row**: *(np.ndarray)*:

    The X matrix row to be populated.

------------
**Returns:**

- tuple: A tuple of arrays:
  - *(np.ndarray)*: updated X_row.
  - *(np.ndarray)*: updated Y_row.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.observation_equations.dz_obs_eq</summary>

---

**dz_obs_eq**(meas_dz: *float*,
              idx_from: *dict*,
              idx_to: *dict*,
              coord_differences: *dict*,
              X_row: *np.ndarray*) -> **Tuple[np.ndarray, np.ndarray]**

---

[[source](src/pysurv/modules/adjustment/observation_equations.py#L179)] Calculates the observation equation for the z-coordinate component of the GNSS vector *(dz)*.

--------------
**Arguments:**

- meas_dz: *(float)*:

    Measured z-coordinate difference.

- **idx_from**: *(dict)*:

    Indices for the station point:

  - **z**: *(int)*: The index of the z-coordinate increment column of the station.

- **idx_to**: *(dict)*:

    Indices for the aim point:

  - **z**: *(int)*: The index of the z-coordinate increment column of the aim point.

- **coord_differences**: *(dict)*:

    Coordinate differences between the points:

  - **dz**: *(float)*: The difference of the z-coordinates of the station and the aim point.

- **X_row**: *(np.ndarray)*:

    The X matrix row to be populated.

------------
**Returns:**

- tuple: A tuple of arrays:
  - *(np.ndarray)*: updated X_row.
  - *(np.ndarray)*: updated Y_row.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.observation_equations.A_obs_eq</summary>

---

**A_obs_eq**(meas_dz: *float*,
              idx_from: *dict*,
              idx_to: *dict*,
              coord_differences: *dict*,
              X_row: *np.ndarray*) -> **Tuple[np.ndarray, np.ndarray]**

---

[[source](src/pysurv/modules/adjustment/observation_equations.py#L215)] Calculates the observation equation for azimuth angle *(A)*.

--------------
**Arguments:**

- meas_A: *(float)*:

    Measured azimuth angle.

- **idx_from**: *(dict)*:

    Indices for the station point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the station.
  - **y**: *(int)*: The index of the y-coordinate increment column of the station.

- **idx_to**: *(dict)*:

    Indices for the aim point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the aim point.
  - **y**: *(int)*: The index of the y-coordinate increment column of the aim point.

- **coord_differences**: *(dict)*:

    Coordinate differences between the points:

  - **dx**: *(float)*: The difference of the x-coordinates of the station and the aim point.
  - **dy**: *(float)*: The difference of the y-coordinates of the station and the aim point.

- **X_row**: *(np.ndarray)*:

    The X matrix row to be populated.

------------
**Returns:**

- tuple: A tuple of arrays:
  - *(np.ndarray)*: updated X_row.
  - *(np.ndarray)*: updated Y_row.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.observation_equations.HZ_obs_eq</summary>

---

**HZ_obs_eq**(meas_dz: *float*,
              idx_from: *dict*,
              idx_to: *dict*,
              coord_differences: *dict*,
              X_row: *np.ndarray*) -> **Tuple[np.ndarray, np.ndarray]**

---

[[source](src/pysurv/modules/adjustment/observation_equations.py#L259)] Calculates the observation equation for horizontal direction *(HZ)*.

--------------
**Arguments:**

- **meas_HZ**: *(float)*:

    Measured horizontal direction.

- **idx_from**: *(dict)*:

    Indices for the station point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the station.
  - **y**: *(int)*: The index of the y-coordinate increment column of the station.
  - **o**: *(int)*: The index of the orientation constant column of the station.

- **idx_to**: *(dict)*:

    Indices for the aim point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the aim point.
  - **y**: *(int)*: The index of the y-coordinate increment column of the aim point.

- **coord_differences**: *(dict)*:

    Coordinate differences between the points:

  - **dx**: *(float)*: The difference of the x-coordinates of the station and the aim point.
  - **dy**: *(float)*: The difference of the y-coordinates of the station and the aim point.
  - **o**: *(float)*: The orientation constant value of the station.

- **X_row**: *(np.ndarray)*:

    The X matrix row to be populated.

------------
**Returns:**

- tuple: A tuple of arrays:
  - *(np.ndarray)*: updated X_row.
  - *(np.ndarray)*: updated Y_row.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.observation_equations.VZ_obs_eq</summary>

---

**VZ_obs_eq**(meas_dz: *float*,
              idx_from: *dict*,
              idx_to: *dict*,
              coord_differences: *dict*,
              X_row: *np.ndarray*) -> **Tuple[np.ndarray, np.ndarray]**

---

[[source](src/pysurv/modules/adjustment/observation_equations.py#L311)] Calculates the observation equation for vertical zenith angle *(VZ)*.

--------------
**Arguments:**

- **meas_VZ**: *(float)*:

    Measured vertical zenith angle.

- **idx_from**: *(dict)*:

    Indices for the station point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the station.
  - **y**: *(int)*: The index of the y-coordinate increment column of the station.
  - **z**: *(int)*: The index of the z-coordinate increment column of the station.

- **idx_to**: *(dict)*:

    Indices for the aim point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the aim point.
  - **y**: *(int)*: The index of the y-coordinate increment column of the aim point.
  - **z**: *(int)*: The index of the z-coordinate increment column of the aim point.

- **coord_differences**: *(dict)*:

    Coordinate differences between the points:

  - **dx**: *(float)*: The difference of the x-coordinates of the station and the aim point.
  - **dy**: *(float)*: The difference of the y-coordinates of the station and the aim point.
  - **dz**: *(float)*: The difference of the z-coordinates of the station and the aim point.

- **X_row**: *(np.ndarray)*:

    The X matrix row to be populated.

------------
**Returns:**

- tuple: A tuple of arrays:
  - *(np.ndarray)*: updated X_row.
  - *(np.ndarray)*: updated Y_row.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.observation_equations.VH_obs_eq</summary>

---

**VH_obs_eq**(meas_dz: *float*,
              idx_from: *dict*,
              idx_to: *dict*,
              coord_differences: *dict*,
              X_row: *np.ndarray*) -> **Tuple[np.ndarray, np.ndarray]**

---

[[source](src/pysurv/modules/adjustment/observation_equations.py#L364)] Calculates the observation equation for vertical horizontal angle *(VH)*.

--------------
**Arguments:**

- **VH_obs_eq**: *(float)*:

    Measured vertical horizontal angle.

- **idx_from**: *(dict)*:

    Indices for the station point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the station.
  - **y**: *(int)*: The index of the y-coordinate increment column of the station.
  - **z**: *(int)*: The index of the z-coordinate increment column of the station.

- **idx_to**: *(dict)*:

    Indices for the aim point:

  - **x**: *(int)*: The index of the x-coordinate increment column of the aim point.
  - **y**: *(int)*: The index of the y-coordinate increment column of the aim point.
  - **z**: *(int)*: The index of the z-coordinate increment column of the aim point.

- **coord_differences**: *(dict)*:

    Coordinate differences between the points:

  - **dx**: *(float)*: The difference of the x-coordinates of the station and the aim point.
  - **dy**: *(float)*: The difference of the y-coordinates of the station and the aim point.
  - **dz**: *(float)*: The difference of the z-coordinates of the station and the aim point.

- **X_row**: *(np.ndarray)*:

    The X matrix row to be populated.

------------
**Returns:**

- tuple: A tuple of arrays:
  - *(np.ndarray)*: updated X_row.
  - *(np.ndarray)*: updated Y_row.

</details>

</details>

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">pysurv.modules.adjustment.robust</summary>

[[source](src/pysurv/modules/adjustment/robust.py)] This module provides various weight functions used in robust estimation.

The module includes four main categories of weight functions:

1. **With Tolerance**: Functions that do not change the weights of observations that fall within a specified range.
2. **Bell Curves**: Functions based on different types of bell-shaped curves.
3. **Trigonometric Functions**: Functions that incorporate trigonometric functions, including hyperbolic functions.
4. **Others**: Miscellaneous functions that do not fall into the previous categories but provide unique methods for adjusting weights.

Each function takes normalized residuals as input and returns an array of reweight factors.
Theoretical values of tuning constants are provided as default.

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Functions</summary>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.huber</summary>

---

**huber**(v: *np.ndarray*, c: *float* = ***1.345***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L48)] Huber M-estimator weight function.

Re-weighting coefficients for residuals are calculated as:
$$1 ~~ for ~~ v \leq c$$
$$\frac{c}{|v|} ~~ for ~~ v \gt c$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(float)*:

    Tuning constant, default is 1.345.

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.slope</summary>

---

**slope**(v: *np.ndarray*, c: *Tuple[float, float]* = ***(2, 2)***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L71)] Slope weight function.

Re-weighting coefficients for residuals are calculated as:
$$1 ~~ for ~~ v \leq c$$
$$1 + \frac{c - |v|}{a} ~~ for ~~ v \gt c and coefficients greater than 0$$
For the rest of them take the value **0**.

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(tuple of two floats)*:

    Tuning constants (c, a), default is (2, 2).

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.hampel</summary>

---

**hampel**(v: *np.array*, c: *Tuple[float, float, float]* = ***(1.7, 3.4, 8.5)***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L97)] Hampel weight function.

Re-weighting coefficients for residuals are calculated as:
$$1 ~~ for ~~ v \leq a$$
$$\frac{a}{|v|} ~~ for ~~ a \lt v \leq b$$
$$\frac{a}{|v|} * \frac{c - |v|}{c - b} ~~ for ~~ b \lt v \leq c$$
$$0 ~~ for ~~ v \gt c$$

--------------
**Arguments:**

- v: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(tuple of three floats)*:

    Tuning constants (a, b, c), default is (1.7, 3.4, 8.5).

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.danish</summary>

---

**danish**(v: *np.ndarray*, c: *float* = ***2.5***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L127)] Danish weight function.

Re-weighting coefficients for residuals are calculated as:
$$1 ~~ for ~~ v \leq c$$
$$exp(-v/c) ~~ for ~~ v \gt c$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(float)*:

    Tuning constant, default is 2.5.

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.epanechnikov</summary>

---

**epanechnikov**(v: *np.ndarray*, c: *Tuple[float, float]* = ***(3.674, 2.0)***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L151)] Epanechnikov weight function.

Re-weighting coefficients for residuals are calculated as:
$$1 - (v/c)^n ~~ for ~~ v \leq c$$
$$0 ~~ for ~~ v \gt c$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(tuple of two floats)*:

    Tuning constants (c, n), default is (3.674, 2.0).

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.tukey</summary>

---

**tukey**(v: *np.ndarray*, c: *Tuple[float, float]* = ***(4.685, 2.0)***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L175)] Tukey weight function.

Re-weighting coefficients for residuals are calculated as:
$$(1 - (v/c)^n)^n ~~ for ~~ v \leq c$$
$$0 ~~ for ~~ v \gt c$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(tuple of two floats)*:

    Tuning constants (c, n), default is (4.685, 2.0).

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.jacobi</summary>

---

**jacobi**(v: *np.ndarray*, c: *Tuple[float, float]* = ***theoretical_c['jacobi']***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L199)] Jacobi weight function.

Re-weighting coefficients for residuals are calculated as:
$$(1 - (v/c)^n)^n * (1 + (v/c)^n)^n ~~ for ~~ v \leq c$$
$$0 ~~ for ~~ v \gt c$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(tuple of two floats)*:

    Tuning constants (c, n), default is (4.687, 1.0).

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.exponential</summary>

---

**exponential**(v: *np.ndarray*, c: *Tuple[float, float]* = ***(2.0, 2.0)***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L223)] Exponential weight function.

Re-weighting coefficients for residuals are calculated as:
$$exp((-v/c)^n) $$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(tuple of two floats)*:

    Tuning constants (c, n), default is (2.0, 2.0).

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.sigma</summary>

---

**sigma**(v: *np.ndarray*, sigma_sq: *float*, c: *Tuple[float, float]* = ***(2.0, 2.0)***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L246)] Sigma weight function.

Re-weighting coefficients for residuals are calculated as:
$$exp((\frac{-v^n}{\sigma_0^2 * c}))$$

The shape of the function changes with each iteration as the value of the residual variance sigma_sq changes.

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **sigma_sq**: *(float)*:

    Residual variance value.

- **c**: *(tuple of two floats)*:

    Tuning constants (c, n), default is (2.0, 2.0).

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.error_func</summary>

---

**error_func**(v: *np.ndarray*, c: *Tuple[float, float]* = ***(1.414, 2.0)***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L272)] Error function weight function.

Re-weighting coefficients for residuals are calculated as:
$$1 - erf((v/c)^n)$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(tuple of two floats)*:

    Tuning constants (c, n), default is (1.414, 2.0).

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.cauchy</summary>

---

**cauchy**(v: *np.ndarray*, c: *Tuple[float, float]* = ***(2.385, 2.0)***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L295)] Cauchy weight function.

Re-weighting coefficients for residuals are calculated as:
$$\frac{1}{1 + (v / c)^n}$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(tuple of two floats)*:

    Tuning constants (c, n), default is (2.385, 2.0).

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.t</summary>

---

**t**(v: *np.ndarray*, k: *int*, c: *Tuple[float, float]* = ***(1.0, 2.0)***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L218)] Student's t weight function.

Re-weighting coefficients for residuals are calculated as:
$$(1 + \frac{v^n}{c * k})^{-\frac{k + 1}{2}}$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **k**: *(int)*:

    Degrees of freedom.

- **c**: *(tuple of two floats)*:

    Tuning constants (c, n), default is (1.0, 2.0).

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.chain_bell</summary>

---

**chain_bell**(v: *np.ndarray*, c: *Tuple[float, float]* = ***(1.0, 1.0)***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L342)] Chain Bell weight function.

Re-weighting coefficients for residuals are calculated as:
$$\frac{1}{cosh((v^n * e)/(2 * c))}$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(tuple of two floats)*:

    Tuning constants (c, n), default is (1.0, 1.0).

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.chain</summary>

---

**chain**(v: *np.ndarray*, c: *float* = ***1.0***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L366)] Chain weight function.

Re-weighting coefficients for residuals are calculated as:
$$-cosh((v * e) / (2 * c)) + 2$$
For coefficients greater than **0**, the rest of them take the value **0**.

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(float)*:

    Tuning constant, default is 1.0.

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.andrews</summary>

---

**andrews**(v: *np.ndarray*, c: *float* = ***4.207***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L391)] Andrews weight function.

Re-weighting coefficients for residuals are calculated as:
$$sinc(v / c) ~~ for ~~ v \leq c$$
$$0 ~~ for ~~ v \gt c$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(float)*:

    Tuning constant, default is 4.207.

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.wave</summary>

---

**wave**(v: *np.ndarray*, c: *float* = ***2.5***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L414)] Wave weight function.

Re-weighting coefficients for residuals are calculated as:
$$\frac{cos((v * \pi )/ c) + 1}{2} ~~ for ~~ v \leq c$$
$$0 ~~ for ~~ v \gt c$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(float)*:

    Tuning constant, default is 2.5.

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.half_wave</summary>

---

**half_wave**(v: *np.ndarray*, c: *float* = ***2.5***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L437)] Half Wave weight function.

Re-weighting coefficients for residuals are calculated as:
$$cos((v * \pi)/(2 * c)) ~~ for ~~ v \leq c$$
$$0 ~~ for ~~ v \gt c$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(float)*:

    Tuning constant, default is 2.5.

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.wigner</summary>

---

**wigner**(v: *np.ndarray*, c: *float* = ***3.137***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L461)] Wigner weight function.

Re-weighting coefficients for residuals are calculated as:
$$\sqrt{1 - (v / c)^2} ~~ for ~~ v \leq c$$
$$0 ~~ for ~~ v \gt c$$

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(float)*:

    Tuning constant, default is 3.137.

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.ellipse_curve</summary>

---

**ellipse_curve**(v: *np.ndarray*, c: *float* = ***2.5***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L484)] Ellipse Curve weight function.

Re-weighting coefficients for residuals are calculated as:
$$\frac{1 + \sqrt{1 - (v/ c)^2}}{2} ~~ for ~~ v \leq c$$
$$\frac{1 - \sqrt{1 - ((v-2*c) / c)^2}}{2} ~~ for ~~ c \gt v \leq 2*c$$
$$0 ~~ for ~~ v \gt 2*c$$

---------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(float)*:

    Tuning constant, default is 2.5.

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">adjustment.robust.trim</summary>

---

**trim**(v: *np.ndarray*, c: *float* = ***2.5***) -> **np.ndarray**

---

[[source](src/pysurv/modules/adjustment/robust.py#L515)] Trim weight function.

Rejects from the set of observations, those whose values of normalized residuals
exceed the value of the tuning constant **c**.

--------------
**Arguments:**

- **v**: *(np.ndarray)*:

    Normalized values of residuals.

- **c**: *(float)*:

    Tuning constant, default is 2.5.

------------
**Returns:**

- **np.ndarray**: Parameters for adjusting the weights.

</details>

</details>

</details>

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.4em; font-weight: bold;">pysurv.modules.Controls</summary>

---

**Controls**(data: *pd.DataFrame*, swap_xy: *bool* = ***False***, *args, **kwargs)

---

[[source](src/pysurv/modules/controls.py#L10)] A specialized **DataFrame subclass** for storing and manipulating control point data, including coordinates (*x*, *y*, *z*) and optionally their standard deviations (*sx*, *sy*, *sz*).

-----------------------------
**Initialization arguments:**

- **data**: *(pd.DataFrame)*:

    The input DataFrame containing control point data.

- **swap_xy**: *(bool)*:

    If True, swaps the *x* and *y* coordinates (and their corresponding sigma values if present).

- **\*args, \*\*kwargs**:

    Additional positional and keyword arguments passed to the pandas DataFrame initializer.

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Properties</summary>

- **coord_labels**: *(pd.Index)*:

    List of only coordinates labels.

- **labels**: *(pd.Index)*:

    List of coordinate labels and orientation constant, if in the set.

- **points_index**: *(dict)*:

    A dictionary mapping point labels (index) to their numerical positions.

- **sigma_labels**: *(pd.Index)*:

    List of only coordinates sigma labels.

- **coord_index**: *(dict)*: 

    A dictionary mapping coordinate labels (index) to their numerical positions.

- **n_unknowns**: *(pd.Index)*: 

    Number of unknown control points.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Methods</summary>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">swap_xy</summary>

---

**swap_xy**()

---

[[source](src/pysurv/modules/controls.py#L118)]  Swap the *x* and *y* coordinates and their corresponding sigma values, if present in the dataset.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">copy_with_type</summary>

---

**copy_with_type**()

---

[[source](src/pysurv/modules/controls.py#L127)] Create and return a **deep copy** of the Controls dataset, preserving the Controls class type.

------------
**Returns:**

- **Controls**: A copy of the current Controls object with its type preserved.

</details>

</details>

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.4em; font-weight: bold;">pysurv.modules.Measurements</summary>

---

**Measurements**(data: *pd.DataFrame*, angle_unit: *str* = ***'grad'***, *args, **kwargs)

---

[[source](src/pysurv/modules/measurements.py#L13)] A specialized pd.DataFrame subclass for storing and manipulating measurements data, and optionally their standard deviations.

-------------------------
**Initialize arguments:**

- **data**: *(pd.DataFrame)*:

    The input DataFrame containing measurements data.

- **angle_unit**: *(str)*:

    Unit of angular measurements in the dataset *('grad', 'gon' or 'deg')*.

- **\*args, \*\*kwargs**:

    Additional positional and keyword arguments passed to the pandas DataFrame initializer.

-----------
**Raises:**

- **ValueError**: If angle_unit is not *'grad'*, *'gon'*, *'deg'* or *'rad'*.

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Properties</summary>

- **types**: (*pd.Index*):

    List of only measurements types labels.

- **sigma**: *(pd.Index)*:

    List of only measurements sigma labels.

- **linear**: *(pd.Index)*:

    List of only linear measurements types labels.

- **linear_sigma**: *(pd.Index)*:

    List of only linear measurements sigma labels.

- **angular**: *(pd.Index)*:

    List of only angular measurements types labels.

- **angular_sigma**: *(pd.Index)*:

    List of only angular measurements sigma labels.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Methods</summary>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">itermeasurements</summary>

---

**itermeasurements**(show_empty: *bool* = ***False***) -> **Tuple[pd.MultiIndex, str, float]**

---

[[source](src/pysurv/modules/measurements.py#L142)] Iterates through all the measurements in the dataset.

--------------
**Arguments:**

- **show_empty**: *(bool, optional)*:

    If True, yields measurements with NaN values. If False, skips NaN values. Default is False.

-----------
**Yields:**

- **tuple**: Measurements in dataset:
  - *(pd.MultiIndex)*: Identifier of points 'FROM' and 'TO'.
  - *(str)*: The measurement type.
  - *(float)*: The measurement value.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">itersigma</summary>

---

**itersigma**(show_empty: *bool* = ***False***) -> **Tuple[pd.MultiIndex, str, float]**

---

[[source](src/pysurv/modules/measurements.py#L167)] Iterates through all the sigmas in the dataset.

--------------
**Arguments:**

- **show_empty**: *(bool, optional)*:

    If True, yields measurements with NaN values. If False, skips NaN values. Default is False.

-----------
**Yields:**

- **tuple**: A tuple containing the index tuple
  - *(pd.MultiIndex)*: Identifier of points 'FROM' and 'TO'.
  - *(str)*: The sigma type.
  - *(float)*: The sigma value.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">iteritems</summary>

---

**iteritems**(show_empty: *bool* = ***False***) -> **Tuple[pd.MultiIndex, str, float]**

---

[[source](src/pysurv/modules/measurements.py#L192)] Iterates through all the measurements and sigma values in the dataset.

--------------
**Arguments:**

- **show_empty**: *(bool, optional)*:

    If True, yields measurements with NaN values. If False, skips NaN values. Default is False.

-----------
**Yields:**

- **tuple**: A tuple containing the index tuple
  - *(pd.MultiIndex)*: Identifier of points 'FROM' and 'TO'.
  - *(str)*: The column type.
  - *(float)*: The value.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">to_disp</summary>

---

**to_disp**(angle_unit: *str* = ***'grad'***)

---

[[source](src/pysurv/modules/measurements.py#L217)] Returns a **deep copy** of the dataset with angles converted to the specified unit.

--------------
**Arguments:**

- **angle_unit**: (*str, optional*):

    The unit in which angles should be displayed. Supported units include 'grad', 'gon' and 'deg'. Default is 'grad'.

------------
**Returns:**

- **Measurements**: A copy of the current dataset with angles converted to the specified unit.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">copy_with_type</summary>

---

**copy_with_type**()

---

[[source](src/pysurv/modules/measurements.py#L239)] Create and return a **deep copy** of the Measurements dataset, preserving the Measurements class type.

------------
**Returns:**

- **Measurements**: A copy of the current Measurements object with its type preserved.

</details>

</details>

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.4em; font-weight: bold;">pysurv.modules.Report</summary>

---

**Report**(results: *dict*, measurements: *Measurements*, path: *str*,
                 methods: *dict* = ***customizations.methods['default']***,
                 report_params: *dict* =  ***customizations.report_params['default']***)

---

[[source](src/pysurv/modules/report/report.py#L25)] A class to generate adjustment report for a surveying control network.

Report contains general information about adjustment computations, sigma evolution plots, control points coordinates,
measurements and standard deviations. Information can be exported to string or to HTML file.

----------------------------------------------------------------------------------------------------------------
**Initialize arguments:**

- **results**: *(dict)*:

    The results from the adjustment process, including matrices, coordinates, residuals, etc.

- **measurements**: *(Measurements)*:

    The measurements dataset used in the adjustment process.

- **path**: *(str)*:

    The path where the report will be saved.

- **methods**: *(dict)*:

    Methods and tuning constants used for adjustment.

- **report_params**: *(dict)*:

    Parameters for customizing the report format.

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Attributes</summary>

- **approx_prec**: *(int)*:

    Precision for the approximate coordinates.

- **adj_prec**: *(int)*:

    Precision for the adjusted coordinates.

- **angle_unit**: *(str)*:

    Unit for angles *('grad', 'gon', or 'deg')*.

- **angle_precision**: *(int)*:

    Precision for angle values.

- **path**: *(str)*:

    The path where the report will be saved.

- **N**: *(np.ndarray)*:

    Matrix of normal equations.

- **cov_b**: *(np.array)*:

    Covariance matrix of the unknowns (adjusted coordinates).

- **cov_Y**: *(np.array)*:

    Covariance matrix of the measurements.

- **cov_r**: *(np.array)*:

    Covariance matrix of the residuals.

- **approx**: *(Controls)*:

    Approximate coordinates.

- **adjusted**: *(Controls)*:

    Adjusted coordinates.

- **obs_residuals**: *(np.array)*:

    Residuals of the observations.

- **obs_residuals_normalized**: *(np.array)*:

    Normalized residuals.

- **obs_sigma_zero**: *(List[float])*:

    List of observation sigma zero values per iteration.

- **pt_sigma_zero**: *(List[float])*:

    List of control points' sigma zero values per iteration.

- **constraints**: *([List[str]], optional)*:

    Applied inner constraints during the adjustment.

- **n_iter**: *(int)*:

    Number of iterations used in the adjustment.

- **n_measurements**: *(int)*:

    Number of measurements in the dataset.

- **n_fixed_coords**: *(int)*:

    Number of fixed reference coordinates.

- **n_sigma_coords**: *(int)*:

    Number of movable reference coordinates.

- **n_unknowns**: *(int)*:

    Number of coordinates and orientation constants to adjust.

- **n_constraints**: *(int)*:

    Number of inner constraints included to the control network.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.2em; font-weight: bold;">Methods</summary>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">to_string</summary>

---

**to_string**(show_plot: *bool* = ***True***) -> **str**

---

[[source](src/pysurv/modules/report/report.py#L329)] Generate a string representation of the report.

--------------
**Arguments**:

- **show_plot**: *(bool, optional)*:

    If True, include the sigma plot in the string. Default is True.

------------
**Returns:**

- **str**: String representation of the report.

</details>

---

<details style="padding-left: 3%;"><summary style="font-size: 1.15em; font-weight: bold;">html</summary>

---

**html**(path: *str or None* = ***None***)

---

[[source](src/pysurv/modules/report/report.py#L352)] Generate an HTML representation of the report and export it to the file.

--------------
**Arguments:**

- **path**: *(str, optional)*:

    The path to save the HTML file. If not provided, default is the path provided at initialization.

</details>

</details>

</details>

</details>

---
