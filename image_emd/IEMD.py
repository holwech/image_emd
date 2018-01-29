import scipy
import numpy as np
import scipy as sp
import scipy.interpolate
import matlab.engine


# Starts the matlab engine and returns the engine as a object
def start_matlab_engine():
    return matlab.engine.start_matlab()


# Performs the Empirical Mode Decomposition on some bidirectional data.
# img is a numpy array of some size MxN and engine is the matlab engine.
# depth is the sensitivity of the extremas. A higher number results in fewer extremas.
def IEMD(img, epsilon, max_imfs, engine, depth=0):
    residue = np.copy(img)
    imfs = np.zeros((max_imfs, img.shape[0], img.shape[1]))
    count = 0
    for i in range(0, max_imfs):
        print("Getting IMF ", count + 1)
        imf = sifting(residue, epsilon, engine, depth)
        residue = residue - imf
        imfs[i,:,:] = imf
        count = count + 1
        print("IMF ", count + 1, " done")
        if monotonic(residue, engine, depth):
            return imfs[:count], residue
    return imfs[:count], residue

# Performs the sifting process until a single IMF is found
def sifting(img, epsilon, engine, depth):
    h_prev = img
    mean = single_sifting(h_prev, engine, depth)
    h_curr = h_prev - mean
    count = 0
    while not (sd(h_curr, h_prev) < epsilon):
        mean = single_sifting(h_curr, engine, depth)
        h_prev = h_curr
        h_curr = h_prev - mean
        count = count + 1
        print("Sifting loop ", count, " sd: ", sd(h_curr, h_prev))
    return h_curr


# Performs a single sifting
def single_sifting(img, engine, depth):
    maxima, minima, maxima_loc, minima_loc = extrema(img, engine, depth)

    x_max, y_max, z_max = triplex_coords(maxima, maxima_loc)
    x_min, y_min, z_min = triplex_coords(minima, minima_loc)

    _,_,_, upper_spline = create_spline(x_max, y_max, z_max)
    _,_,_, lower_spline = create_spline(x_min, y_min, z_min)

    upper_img = reconstruct_image(img.shape, upper_spline)
    lower_img = reconstruct_image(img.shape, lower_spline)

    mean = surface_mean(upper_img, lower_img)

    return mean

# Checks whether the residue is monotonic or not, that is, if the residue has less than 2 extrema.
def monotonic(residue, engine, depth):
    _,_, maxima_loc, minima_loc = extrema(residue, engine, depth)
    num_extrema = np.sum(maxima_loc) + np.sum(minima_loc)
    if num_extrema < 2:
        return True
    return False


# Calculates the stopping criterion value (sd = standard deviation??)
def sd(h_curr, h_prev):
    numerator = np.square(h_prev - h_curr)
    denominator = np.square(h_prev)
    return np.sum(numerator / denominator)

# Calculated the piecewise mean of the upper and lower interpolation
def surface_mean(upper_img, lower_img):
    return (upper_img + lower_img) / 2

# Reconstructs the upper and lower plate based on found spline.
def reconstruct_image(size, spline):
    new_img = np.zeros(size)
    for iy, row in enumerate(new_img):
        for ix, col in enumerate(row):
            new_img[iy, ix] = spline(ix, iy)
    return new_img

# Converts a 2D array to its (x, y, z) coordinates in 3D space
def triplex_coords(locations, values):
    x = np.array([])
    y = np.array([])
    z = np.array([])
    for iy, row in enumerate(locations):
        for ix, col in enumerate(row):
            if col:
                y = np.append(y, iy)
                x = np.append(x, ix)
                z = np.append(z, values[iy, ix])
    return x, y, z

# Returns the values of the mixima and minima in their respective locations
def extrema(img, engine, depth):
    mat_img = matlab.double(img.tolist())
    maxima_loc = mat2np(engine.imextendedmax(mat_img, depth))
    minima_loc = mat2np(engine.imextendedmin(mat_img, depth))
    maxima = maxima_loc * img
    minima = minima_loc * img
    return maxima, minima, maxima_loc, minima_loc

# A fast conversion from matlab.double array to numpy array
def mat2np(mat):
    return np.array(mat._data).reshape(mat.size, order='F')

def create_spline(x,y,z):
    x_grid = np.linspace(0, 75, len(x))
    y_grid = np.linspace(0, 75, len(y))
    X, Y = np.meshgrid(x_grid, y_grid, indexing='xy')
    Z = np.zeros((x.size, z.size))
    spline = sp.interpolate.Rbf(x,y,z,function='thin-plate')
    Z = spline(X,Y)
    return X, Y, Z, spline
