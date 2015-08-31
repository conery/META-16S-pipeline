# Paths to applications and resource files, other configuration info

# The resources directory contains the reference sequences used by 
# the chimera filtering script and map-to-refernce script

import os
import sys

resource_dir = os.path.join(sys.path[0], '..', 'resources')

# Define the path to the jar file with the executable of the RDP Classifier

classifier_path = '~/Applications/rdp_classifier_2.10.2/dist/classifier.jar'

