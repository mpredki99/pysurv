# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

import re

# ======================================================================================
# Regular expression pattern to match comment lines in data files.
# Lines starting with '#' are considered comments.
# ======================================================================================
COMMENT = re.compile(r"^[\s\t]*#")
