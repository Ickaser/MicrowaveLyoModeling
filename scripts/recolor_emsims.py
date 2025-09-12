import numpy as np
import scipy
import os
import matplotlib
import copy
import matplotlib.pyplot as plt
from matplotlib.pyplot import imread, imsave

def_colorbar_loc = [[909, 22], [953 - 20, 959]]
def_working_loc = [[95, 189], [857, 708]]


def recolor_image(
    filename='image.png',
    colorbar_loc=def_colorbar_loc,
    working_loc=def_working_loc,
    colorbar_orientation='auto',
    colorbar_direction=-1,
    new_cmap='viridis',
    normalize_before_compare=False,
    max_rgb='auto',
    threshold=0.4,
    saturation_threshold=0.25,
    compare_hue=True,
    show_plot=True,
    debug=False,
):

    """
    This script reads in an image file (like .png), reads the image's color bar (you have to tell it where), interprets
    the color map used in the image to convert colors to values, then recolors those values with a new color map and
    regenerates the figure. Useful for fixing figures that were made with rainbow color maps.

    Parameters

    -----------

    :param filename: Full path and filename of the image file.
    :param colorbar_loc: Location of color bar, which will be used to analyze the image and convert colors into values.
        Pixels w/ 0,0 at top left corner: [[left, top], [right, bottom]]
    :param working_loc: Location of the area to recolor. You don't have to recolor the whole image.
        Pixels w/ 0,0 at top left corner: [[left, top], [right, bottom]], set to [[0, 0], [-1, -1]] to do everything.
    :param colorbar_orientation: Set to 'x', 'y', or 'auto' to specify whether color map is horizontal, vertical,
        or should be determined based on the dimensions of the colorbar_loc
    :param colorbar_direction: Controls direction of ascending value
        +1: colorbar goes from top to bottom or left to right.
        -1: colorbar goes from bottom to top or right to left.
    :param new_cmap: String describing the new color map to use in the recolored image.
    :param normalize_before_compare: Divide r, g, and b each by (r+g+b) before comparing.
    :param max_rgb: Do the values of r, g, and b range from 0 to 1 or from 0 to 255? Set to 1, 255, or 'auto'.
    :param threshold: Sum of absolute differences in r, g, b values must be less than threshold to be valid
        (0 = perfect, 3 = impossibly bad). Higher numbers = less chance of missing pixels but more chance of recoloring
        plot axes, etc.
    :param saturation_threshold: Minimum color saturation below which no replacement will take place
    :param compare_hue: Use differences in HSV instead of RGB to determine with which index each pixel should be
        associated.
    :param show_plot: T/F: Open a plot to explain what is going on. Also helpful for checking your aim on the colorbar
        coordinates and debugging.
    :param debug: T/F: Print debugging information.
    """

    def printd(string_in):
        """
        Prints debugging statements
        :param string_in: String to print only if debug is on.
        :return: None
        """
        if debug:
            print(string_in)
        return

    print('Recoloring image: {:} ...'.format(filename))

    # Determine tag name and load original file into the tree
    fn1 = filename.split(os.sep)[-1]  # Filename without path
    fn2 = fn1.split(os.extsep)[0]  # Filename without extension (so new filename can be built later)
    ext = fn1.split(os.extsep)[-1]  # File extension
    path = os.sep.join(filename.split(os.sep)[0:-1])  # Path; used later to save results.
    a = imread(filename).astype(np.float32)
    printd(f'Read image; shape = {np.shape(a)}')

    if max_rgb == 'auto':
        # Determine if values of R, G, and B range from 0 to 1 or from 0 to 255
        if a.max() > 1:
            max_rgb = 255.0
        else:
            max_rgb = 1.0
    # Normalize a so RGB values go from 0 to 1 and are floats.
    a /= max_rgb

    # Extract the colorbar
    x = np.array([colorbar_loc[0][0], colorbar_loc[1][0]])
    y = np.array([colorbar_loc[0][1], colorbar_loc[1][1]])
    cb = a[y[0]:y[1], x[0]:x[1]]
    # Take just the working area, not the whole image
    xw = np.array([working_loc[0][0], working_loc[1][0]])
    yw = np.array([working_loc[0][1], working_loc[1][1]])
    a1 = a[yw[0]:yw[1], xw[0]:xw[1]]

    # Pick color bar orientation
    if colorbar_orientation == 'auto':
        if np.diff(x) > np.diff(y):
            colorbar_orientation = 'x'
        else:
            colorbar_orientation = 'y'
        printd('Auto selected colorbar_orientation')
    printd('Colorbar orientation is {:}'.format(colorbar_orientation))

    # Analyze the colorbar
    if colorbar_orientation == 'y':
        cb = np.nanmean(cb, axis=1)
    else:
        cb = np.nanmean(cb, axis=0)
    if colorbar_direction < 0:
        cb = cb[::-1]
    # Compress colorbar to only count unique colors
    # If the array gets too big, it will fill memory and crash python: https://github.com/numpy/numpy/issues/14136
    dcb = np.append(1, np.sum(abs(np.diff(cb[:, 0:3], axis=0)), axis=1))
    cb = cb[dcb > 0]

    # ISW edit: trim colorbar to a maximum of 128 unique colors
    if cb.shape[0] > 128:
        # orig_num = cb.shape[0]
        # step = orig_num // 128
        # np.array(
        # cols = [0] + cb.shape
        cb = np.concatenate((cb[0::(cb.shape[0] // 128), :], cb[-1:, :]), axis=0)
    printd('Compressed colorbar to {:} unique colors'.format(np.shape(cb)[0]))

    # Find and mask of special colors that should not be recolored
    n1a = np.sum(a1[:, :, 0:3], axis=2)
    replacement_mask = np.ones(np.shape(n1a), bool)
    for col in [0, 3]:  # Black and white will come out as 0 and 3.
        mask_update = n1a != col
        if mask_update.max() == 0:
            print('Warning: masking to protect special colors prevented all changes to the image!')
        else:
            printd('Good: Special color mask {:} allowed at least some changes'.format(col))
        replacement_mask *= mask_update
        if replacement_mask.max() == 0:
            print('Warning: replacement mask will prevent all changes to the image! '
                  '(Reached this point during special color protection)')
        printd('Sum(replacement_mask) = {:}    (after considering special color {:})'
               .format(np.sum(np.atleast_1d(replacement_mask)), col))
    # Also apply limits to total r+g+b
    # replacement_mask *= n1a > 0.25
    # replacement_mask *= n1a < 2.75
    if replacement_mask.max() == 0:
        print('Warning: replacement mask will prevent all changes to the image! '
              '(Reached this point during total r+g+b+ limits)')
    printd('Sum(replacement_mask) = {:}    (after considering r+g+b upper threshold)'
           .format(np.sum(np.atleast_1d(replacement_mask))))
    if saturation_threshold > 0:
        hsv1 = matplotlib.colors.rgb_to_hsv(a1[:, :, 0:3])
        sat = hsv1[:, :, 1]
        printd('Saturation ranges from {:} <= sat <= {:}'.format(sat.min(), sat.max()))
        sat_mask = sat > saturation_threshold
        if sat_mask.max() == 0:
            print('Warning: saturation mask will prevent all changes to the image!')
        else:
            printd('Good: Saturation mask will allow at least some changes')
        replacement_mask *= sat_mask
        if replacement_mask.max() == 0:
            print('Warning: replacement mask will prevent all changes to the image! '
                  '(Reached this point during saturation threshold)')

    printd(f'shape(a1) = {np.shape(a)}')
    printd(f'shape(cb) = {np.shape(cb)}')
    # Find where on the colorbar each pixel sits
    if compare_hue:
        # Difference in hue
        hsv1 = matplotlib.colors.rgb_to_hsv(a1[:, :, 0:3])
        hsv_cb = matplotlib.colors.rgb_to_hsv(cb[:, 0:3])
        d2 = abs(hsv1[:, :, :, np.newaxis] - hsv_cb.T[np.newaxis, np.newaxis, :, :])
        # d2 = d2[:, :, 0, :]  # Take hue only
        d2 = np.sum(d2, axis=2)
        printd('  shape(d2) = {:}    (hue version)'.format(np.shape(d2)))
    else:
        # Difference in RGB
        if normalize_before_compare:
            # Difference of normalized RGB arrays
            n1 = n1a[:, :, np.newaxis]
            n2 = np.sum(cb[:, 0:3], axis=1)[:, np.newaxis]
            w1 = n1 == 0
            w2 = n2 == 0
            n1[w1] = 1
            n2[w2] = 1
            d = (a1/n1)[:, :, 0:3, np.newaxis] - (cb/n2).T[np.newaxis, np.newaxis, 0:3, :]
        else:
            # Difference of non-normalized RGB arrays
            d = (a1[:, :, 0:3, np.newaxis] - cb.T[np.newaxis, np.newaxis, 0:3, :])

        printd(f'Shape(d) = {np.shape(d)}')

        d2 = np.sum(np.abs(d[:, :, 0:3, :]), axis=2)  # 0:3 excludes the alpha channel from this calculation
    printd('Processed colorbar')

    index = d2.argmin(axis=2)
    md2 = d2.min(axis=2)
    index_valid = md2 < threshold
    if index_valid.max() == 0:
        print('Warning: minimum difference is greater than threshold: all changes rejected!')
    else:
        printd('Good: Minimum difference filter is lower than threshold for at least one pixel.')
    printd('Sum(index_valid) = {:}     (before *= replacement_mask)'.format(np.sum(np.atleast_1d(index_valid))))
    printd('Sum(replacement_mask) = {:}    (final, before combining w/ index_valid)'
           .format(np.sum(np.atleast_1d(replacement_mask))))
    index_valid *= replacement_mask
    if index_valid.max() == 0:
        print('Warning: index_valid mask prevents all changes to the image after combination w/ replacement_mask.')
    else:
        printd('Good: Mask will allow at least one pixel to change.')
    printd('Sum(index_valid) = {:}'.format(np.sum(np.atleast_1d(index_valid))))
    value = index/(len(cb)-1.0)
    printd('Index ranges from {:} to {:}'.format(index.min(), index.max()))

    # Make a new image with replaced colors
    b = matplotlib.cm.ScalarMappable(cmap=new_cmap).to_rgba(value)  # Remap everything
    printd('shape(b) = {:}, min(b) = {:}, max(b) = {:}'.format(np.shape(b), b.min(), b.max()))
    c = copy.copy(a1)  # Copy original
    c[index_valid, 0:3] = b[index_valid, 0:3]  # Transfer only pixels where color was close to colormap

    # Transfer working area to full image
    c2 = copy.copy(a)  # Copy original full image
    c2[yw[0]:yw[1], xw[0]:xw[1], :] = c  # Replace working area
    # c2[:, :, 3] = a[:, :, 3]  # Preserve original alpha channel

    # Save the image in the same path as the original but with _recolored added to the filename.
    new_filename = '{:}{:}{:}_recolored{:}{:}'.format(path, os.sep, fn2, os.extsep, ext)
    imsave(new_filename, c2)

    print('Done recoloring. Result saved to {:} .'.format(new_filename))

    if show_plot:
        # Setup figure for showing things to the user
        f, axs = plt.subplots(2, 3)
        axo = axs[0, 0]  # Axes for original figure
        axoc = axs[0, 1]  # Axes for original color bar
        axf = axs[0, 2]  # Axes for final figure
        axm = axs[1, 1]  # Axes for mask
        axre = axs[1, 2]  # Axes for recolored section only (it might not be the whole figure)
        axraw = axs[1, 0]  # Axes for raw recoloring result before masking

        for ax in axs.flatten():
            ax.set_xlabel('x pixel')
            ax.set_ylabel('y pixel')

        axo.set_title('Original image w/ colorbar ID overlay')
        axoc.set_title('Color progression from original colorbar')
        axm.set_title('Mask')
        axre.set_title('Recolored section')
        axraw.set_title('Raw recolor result (no masking)')
        axf.set_title('Final image')

        axoc.set_xlabel('Index')
        axoc.set_ylabel('Value')

        # Show the user where they placed the color bar and working location
        axo.imshow(a)
        xx = x[np.array([0, 0, 1, 1, 0])]
        yy = y[np.array([0, 1, 1, 0, 0])]
        axo.plot(xx, yy, '+-', label='colorbar')

        xxw = xw[np.array([0, 0, 1, 1, 0])]
        yyw = yw[np.array([0, 1, 1, 0, 0])]
        axo.plot(xxw, yyw, '+-', label='target')

        tots = np.sum(cb[:, 0:3], axis=1)

        if normalize_before_compare:
            # Normalized version
            axoc.plot(cb[:, 0] / tots, 'r', label='r/(r+g+b)', lw=2)
            axoc.plot(cb[:, 1] / tots, 'g', label='g/(r+g+b)', lw=2)
            axoc.plot(cb[:, 2] / tots, 'b', label='b/(r+g+b)', lw=2)
            axoc.set_ylabel('Normalized value')
        else:
            axoc.plot(cb[:, 0], 'r', label='r', lw=2)
            axoc.plot(cb[:, 1], 'g', label='g', lw=2)
            axoc.plot(cb[:, 2], 'b', label='b', lw=2)
            # axoc.plot(cb[:, 3], color='gray', linestyle='--', label='$\\alpha$')
            axoc.plot(tots, 'k', label='r+g+b')

        # Display the new colors with no mask, the mask, and the recolored section
        axraw.imshow(b)
        axm.imshow(index_valid)
        axre.imshow(c)

        # Display the final result
        axf.imshow(c2)

        # Finishing touches on plots
        axo.legend(loc=0).set_draggable(True)
        axoc.legend(loc=0).set_draggable(True)
        plt.show()

    return

cmap = "plasma"
thresh = 0.4
sat_thresh = 0.4

colorbar_loc = [[2250, 494], [2300, 2784]]
sec1 = [[380, 450], [1900, 1400]]
sec2 = [[380, 1950], [1900, 2900]]
sec3 = [[2200, 475], [2330, 2880]]

recolor_image(os.path.abspath("../data/emsims/glass_all.png"), colorbar_loc, sec1,
              "y", -1, cmap, False, "auto", thresh, sat_thresh, False, False, True)
recolor_image(os.path.abspath("../data/emsims/glass_all_recolored.png"), colorbar_loc, sec2,
              "y", -1, cmap, False, "auto", thresh, sat_thresh, False, False, True)
recolor_image(os.path.abspath("../data/emsims/glass_all_recolored_recolored.png"), colorbar_loc, sec3,
              "y", -1, cmap, False, "auto", thresh, sat_thresh, False, False, True)

colorbar_loc = [[2254, 448], [2283, 2514]]
sec1 = [[275, 890], [1975, 1280]]
sec2 = [[260, 2200], [1940, 2500]]
sec3 = [[2220, 430], [2320, 2550]]


recolor_image(os.path.abspath("../data/emsims/ice_all.png"), colorbar_loc, sec1,
              "y", -1, cmap, False, "auto", thresh, sat_thresh, False, False, True)
recolor_image(os.path.abspath("../data/emsims/ice_all_recolored.png"), colorbar_loc, sec2,
              "y", -1, cmap, False, "auto", thresh, sat_thresh, False, False, True)
recolor_image(os.path.abspath("../data/emsims/ice_all_recolored_recolored.png"), colorbar_loc, sec3,
              "y", -1, cmap, False, "auto", thresh, sat_thresh, False, False, True)