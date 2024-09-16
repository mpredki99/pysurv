# Coding: UTF-8

# Copyright (C) 2024 Michał Prędki
# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files in the repository.

import base64
from io import BytesIO
import matplotlib.pyplot as plt
import pandas as pd


def to_html(plot: plt.figure, general_info: pd.DataFrame,
            controls_table: pd.DataFrame, measurements_table: pd.DataFrame, footer: list) -> str:
    """
    Built HTML representation of report.

    --------------------------------------------------------------------------------------------------------------------
    Arguments:

    - plot: (figure): Plot of sigma zero values.
    - general_info: (pd.DataFrame): DataFrame with general information about adjustment.
    - controls_table: (pd.DataFrame): DataFrame containing information about control points.
    - measurements_table: (pd.DataFrame): DataFrame containing information about measurements.
    - footer: (list): Footer of the report.
    """

    # Convert plot to image
    file = BytesIO()
    plot.savefig(file, format='png')
    encoded_image = base64.b64encode(file.getvalue()).decode('utf-8')
    image_tag = f"<img src='data:image/png;base64,{encoded_image}' class='plot-image'>"
    # CSS and JavaScript content
    bootstrap_css = """
    <!-- Bootstrap CSS -->
    <link href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css" rel="stylesheet">
    """
    dark_mode_toggle = """
    <button id="darkModeToggle">Toggle Dark Mode</button>
    <script>
        const toggle = document.getElementById('darkModeToggle');
        toggle.addEventListener('click', () => {
        document.body.classList.toggle('dark-mode');
    });
    </script>
    <style>
        .dark-mode {
            background-color: #2c2c2c;
            color: #f1f1f1;
        }
        .dark-mode table {
            background-color: #3e3e3e;
            color: #f1f1f1;
        }
        .dark-mode .plot-image {
            filter: invert(1) hue-rotate(180deg);
        }
    </style>
    """
    print_mode_styles = """
    <style>
    @media print {
            body {
                    background: none;
                    color: black;
                }
            table {
                border: 1px solid black;
            }
        }
    </style>
    """
    report_styles = """
    <style>
        html {
            scroll-behavior: smooth;
        }
        body {
            font-family: 'Cambria';
        }
        h3, h5 {
            font-weight: 700;
        }
        p {
            margin-left: 5%;
        }
        /* Tables styles */
        table {
            border: 2px solid #000000;
            border-collapse: collapse;
        }
        table th, table td {
            border: 1px solid #000000;
            background-color: rgba(0,150,30,.1);
        }
        .half-width-table {
            width: 50%;
            margin: 0 auto;
        }
        .full-width-table {
            width: 90%;
            margin: 0;
        }
        .full-width-table table th,
        .full-width-table table td {
            text-align: center;
        }
        .full-width-table table th {
            border-bottom: 2px solid #000000;
            border-top: 2px solid #000000;
            background-color: rgba(0,150,30,.5);
        }
        .full-width-table:nth-of-type(2) table td:nth-child(1),
        .full-width-table:nth-of-type(3) table td:nth-child(1),
        .full-width-table:nth-of-type(3) table td:nth-child(2) {
            font-weight: bold;
        }
        .table-striped tbody tr:nth-of-type(odd) {
            background-color: rgba(0,150,30,.2);
        }
        .table-hover tbody tr:hover {
            background-color: rgba(0,150,30,.4);
        }
    </style>
    """
    # CSS classes for to_html method
    table_style = 'table table-striped table-hover table-sm table responsive'
    # Generate HTML content
    html_content = f"""
<!DOCTYPE html>
    <html>
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            {bootstrap_css}
            {dark_mode_toggle}
            {print_mode_styles}
            {report_styles}
            <title>Least Squares Adjustment Results</title>
        </head>
        <body>
            <center>
                <h3>Least Squares Adjustment Results</h3>

                <h5>General Information</h5>
                <div class="half-width-table">
                    {general_info.to_html(header=False, na_rep='-', classes=table_style)
        .replace('obs_sigma_zero', '&sigma;<sub>0</sub>')
        .replace('pt_sigma_zero', '&sigma;<sub>&beta;</sub>')
        .replace('NaN', '')
        }
            </div>

            {image_tag}

            <h5>Control Points</h5>
            <div class="full-width-table">
                {controls_table.reset_index().to_html(index=False, na_rep='', classes=table_style)
        .replace('_0_', '<sub>0</sub>').replace('_s_', '&sigma;').replace('NaN', '')
        }
            </div>

            <h5>Measurements</h5>
            <div class="full-width-table">
                {measurements_table.reset_index().to_html(index=False, na_rep='', classes=table_style)
        .replace('_0_', '<sub>0</sub>').replace('_s_', '&sigma;').replace('NaN', '')
        }
            </div>
        </center>
        <p>{'</p><p>'.join(footer)}
    </body>
</html>
        """

    return html_content
