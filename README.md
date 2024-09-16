# pysurv

<div style="text-align: justify">

**pysurv** is a Python package for adjusting surveying control networks. The package supports importing data from CSV files and performing ordinary, weighted, and robust least squares adjustment. It also allows the free adjustment approach to be combined with ordinary, weighted, and robust methods. Additionally, you can mix these methods when adjusting observations and reference points. After completing the calculations, a HTML report with the adjustment results can be generated.

## Features

- Import surveying data from standardized CSV files.
- Perform ordinary, weighted or robust least squares adjustment of control network.
- Generate an report with adjustment results.

## Installation

You can install the package via `pip` directly from Git-Hub repository:

``` bash
pip install git+https://github.com/mpredki99/pysurv.git
```

## Usage

### Project Setup

1. **Create Project Folder:** Create a project folder on your computer where you will place the files containing the approximate coordinates of control points and measurement data.
2. **Prepare Input Data Files:** Place the following CSV files in the project folder:

    2.1 **controls.csv**: This file should contain the approximate coordinates of the control points for adjustment. The columns in the file should follow this naming convention:
    - ***ID***: The unique identifier of the point.
    - ***x***: The x-coordinate of control point.
    - ***y***: The y-coordinate of control point.
    - ***z***: The z-coordinate of control point.

    The file may also include optional columns for the estimated **standard deviations** of the control point coordinates: ***sx***, ***sy***, ***sz***. If the standard deviation values are unknown, leave these cells blank. The program will use a default value to determine the point's weight in that case. Setting a value of *-1* in these columns will assign a weight of the point coordinate to *0*, in practice, excluding the point as a reference.

    Alternatively, you can use *easting* (***E***), *northing* (***N***), and *elevation* (***el***) or *height* (***H***) coordinates. The program will automatically convert these into ***x***, ***y***, and ***z*** coordinates during import.

    **Note**: By default, the program maps *easting* (***E***) to ***y*** and *northing* (***N***) to ***x***.

    <div align=center>

    **Example of *"controls.csv"***

    | ID  |   x     |    y    |    z   |
    |:---:|:-------:|:-------:|:------:|
    | P01 | 317.08  | 1606.89 | 103.57 |
    | P02 | 706.43  | 1785.03 | 99.73  |
    | P03 | 1063.29 | 1552.50 | 99.19  |
    | P04 | 852.32  | 1293.80 | 97.48  |

    </div>

    2.2 **measurements.csv**: This file should contain the surveying measurement data with the following columns:
    - ***FROM***: The identifier of the point **from** which the measurement was taken.
    - ***TO***: The identifier of the point **to** which the measurement was taken.
    - ***SD***: The 3D spatial distance.
    - ***HD***: The 2D horizontal distance.
    - ***VD***: The 1D vertical distance.
    - ***dx***: The x-component of the GNSS vector in the **plane** coordinate system.
    - ***dy***: The y-component of the GNSS vector in the **plane** coordinate system.
    - ***dz***: The z-component of the GNSS vector in the **plane** coordinate system.
    - ***A***: The azimuth angle.
    - ***HZ***: The horizontal direction.
    - ***VZ***: The vertical zenith angle.
    - ***VH***: The vertical horizontal angle.

    You can also include estimated **standard deviations** for the measurements by creating corresponding columns with the prefix *"s"* (e.g., ***sSD***, ***sHD***, ***sHZ***, ***sVZ***). If the standard deviation values are unknown, leave these cells blank, and the program will use default values to determine the observation weight.

    <div align=center>

    **Example of *"measurements.csv"***

    | FROM  |  TO   |   HD    |  sHD   |    HZ    |    VZ    |
    |:-----:|:-----:|:-------:|:------:|:--------:|:--------:|
    |  P01  |  P02  | 428.173 | 0.0024 | 308.7758 | 100.5715 |
    |  P01  |  P03  | 748.195 | 0.0027 | 276.8270 | 100.3728 |
    |  P02  |  P01  |         |        | 384.1063 | 99.4291  |
    |  P02  |  P03  | 425.933 | 0.0024 | 120.0253 | 100.0820 |
    |  P02  |  P04  | 512.430 | 0.0025 | 75.1682  | 100.2799 |
    |  P03  |  P01  |         |        | 27.2014  | 99.6268  |
    |  P03  |  P04  | 333.814 | 0.0023 | 88.2790  | 100.3253 |
    |  P04  |  P01  | 620.088 | 0.0026 | 270.0645 | 99.3742  |
    |  P04  |  P02  |         |        | 222.1380 | 99.7206  |
    |  P04  |  P03  |         |        | 160.2063 | 99.6756  |

    </div>

### Performing adjustment in your project

After installing the package, you can create an instance of the Project class and run the adjustment as follows:

``` python
>>> from pysurv import Project

>>> # Create a project instance
>>> project = Project('path/to/project_folder')

>>> # Perform the adjustment and export an HTML report
>>> project.adjust(report_format='html')

```

#### Project Customisations

1. **Methods**

    In the example above, the adjustment is performed using the default parameters, i.e., the weighted least squares method for adjusting observations without free adjustment approach involve. To view the current default methods values, you can retrieve them from the customisations module.

    ``` python
    >>> from pysurv import customizations
    
    >>> print(customizations.methods)
   
    ```

    Output:

    ```output
    >>> {'default': {'observations': 'weighted', 
    ...              'obs_C': None, 
    ...              'free': None, 
    ...              'free_C': None}
    ...  }
    ```

    The code above shows how to check the contents of the **methods** dictionary. To add custom parameters to the dictionary, use the **update()** method.

    ```python
    >>> from pysurv import customizations

    >>> # Add custom parameters to the methods dictionary
    >>> customizations.methods.update({'user': {'observations': 'huber',
    ...                                         'obs_C': 1.55,
    ...                                         'free': 'tukey',
    ...                                         'free_C': None}})
    >>> # Display the updated dictionary
    >>> print(customizations.methods)
   
    ```

    Output:

    ```output
    >>> {'default': {'observations': 'weighted', 
    ...             'obs_C': None, 
    ...              'free': None, 
    ...              'free_C': None}, 
    ...  'user': {'observations': 'huber', 
    ...           'obs_C': 1.55, 
    ...           'free': 'tukey', 
    ...           'free_C': None}
    ...  }
    ```

    In the example above, the custom user-defined method specifies the use of Huber's robust M-estimator to adjust the measurements, with a tuning constant of 1.55. The free adjustment will be performed using an inner constraint matrix, and the point weights will be adjusted using robust Tukey's M-estimator. Since the tuning constant for the free adjustment is not defined by the user (**None**), the default theoretical value for this method will be used.

    After adding the custom parameters to the dictionary, the adjustment can be performed according to the user-defined parameters.

    ```python
    >>> from pysurv import Project
    
    >>> # Create a project instance with the user's defined method and tuning constants
    >>> project = Project('path/to/project_folder', methods='user')
    
    >>> # Perform the adjustment and export an HTML report
    >>> project.adjust(report='html')
    
    ```

2. **Controls and measurements standard deviations**

    In addition to defining the methods and tuning constants for performing the adjustment in the project, you can also customize the default values for the **standard deviations** of the **measurements** and **control point locations**, which will be used to determine the default weights if they are not specified in the input files. These values can be accessed from the customisations module.

    ```python
    >>> from pysurv import customizations

    >>> print(customizations.default_measurement_sigma)
   
    ```

    Output:

    ```output
    >>> {'default': {'sSD': (0.003, 0.002),
    ...              'sHD': (0.003, 0.002),
    ...              'sVD': (0.003, 0.002),
    ...              'sdx': (0.002, 0.001),
    ...              'sdy': (0.002, 0.001),
    ...              'sdz': (0.002, 0.001),
    ...              'sHZ': 0.0020,
    ...              'sVZ': 0.0020,
    ...              'sVH': 0.0020,
    ...              'sA': 0.0020,
    ...              'sP': 0.05}
    ... }
    ```

    The standard deviations of linear measurements are presented in two-element tuples because they contain a **constant** error component and a **distance-dependent** error (*p.p.m.*).

3. **Iteration process parameters**

    The next parameters that can be customized are the iteration process parameters. These parameters define:

    - **tolerance**: the threshold at which further adjustments to the control point coordinates become negligible, allowing the calculations to be terminated
    - **max_iter**: maximum number of iterations to be performed during calculations.

    ```python
    >>> from pysurv import customizations

    >>> print(customizations.iterate_params)
   
    ```

    Output:

    ```output
    >>> {'default': {'tolerance': 0.0001, 
    ...              'max_iter': 100}
    ... }
    ```

4. **Report parameters**

    The final set of parameters allows you to customize the appearance of the report. These parameters determine the number of decimal places used for displaying approximate and adjusted values of control point coordinates and linear measurements. Additionally, they define the units in which angles are displayed and the number of decimal places for angular measurements.

    ```python
    >>> from pysurv import customizations
    
    >>> print(customizations.iterate_params)
        
    ```

    Output:

    ```output
    >>> {'default': {'approx_precision': 3, 
    ...              'adjusted_precision': 4, 
    ...              'angle_unit': 'grad', 
    ...              'angle_precision': 4}
    ... }
    ```

5. **Running customized project**

    After defining all the parameters of the adjustment and report, you can create an instance of the Project class including user-defined parameters.

    ```python
    >>> from pysurv import Project
    
    >>> # Create a project instance with the user's defined method and tuning constants
    >>> project = Project('path/to/project_folder', 
    ...                    methods='user',
    ...                    measurement_errors='user',
    ...                    iterate_params='user', 
    ...                    report_params='user')
    
    >>> # Perform the adjustment and export an HTML report
    >>> project.adjust(report='html')
    
    ```

    Additional parameters possible to set when creating a project instance:
    - **angle_unit** - specifies the unit in which the angular measurements in the measurement data set (measurements.csv) are given. This parameter takes the values *'grad'* or *'gon'* for measurements expressed in **gradians**, and *'deg'* for measurements expressed in **decimal degrees**. Defaults to *'grad'*.
    - **swap_xy** - determines whether to swap *x* and *y* coordinate values when importing control point coordinates. Dafaults to *False*.

    ```python
    >>> from pysurv import Project
    
    >>> # Create a project instance with the user's defined method and tuning constants
    >>> project = Project('path/to/project_folder', 
    ...                    methods='user',
    ...                    measurement_errors='user',
    ...                    swap_xy=False, 
    ...                    angle_unit='grad' 
    ...                    iterate_params='user', 
    ...                    report_params='user')
    
    >>> # Perform the adjustment and export an HTML report
    >>> project.adjust(report='html')
    
    ```

### Using pysurv functionality in your project

The Project class allows you to easily import data and perform surveying control network adjustment, but the entire ***pysurv*** package contains a lot of functionalities that you may find useful during the development of your projects.

To import them into your project simply, use the import statement.

```python
import pysurv.modules as ps
```

For more detailed documentation of pysurv functionalities see [doc.md](doc.md).

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

</div>
