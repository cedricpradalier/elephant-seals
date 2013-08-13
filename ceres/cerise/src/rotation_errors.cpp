
#include "rotation_errors.h"

double cerise::AccelerometerError::G[3] = {0,0,-10};

// Extracted with igrf, from the matlab side, for the beach

double cerise::MagnetometerError::B[3] = {9.6920e+03, -1.4563e+04, -4.5355e+04};
