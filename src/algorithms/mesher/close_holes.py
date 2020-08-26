import numpy as np
import matplotlib.pyplot as plt
import skimage.io
import skimage.morphology, skimage.data

def close_holes(img, flow_pixel, no_flow_pixel):
    labels = skimage.morphology.label(img, background = no_flow_pixel)

    unique_in, i_in = np.unique(labels[:, 0], return_index=True)
    unique_out, i_out = np.unique(labels[:, -1], return_index=True)

    img_pixel_in = np.where(img[i_in, 0] == flow_pixel)[0]
    img_pixel_out = np.where(img[i_out, -1] == flow_pixel)[0]

    print(unique_in[img_pixel_in])
    print( unique_out[img_pixel_out])
    print(unique_in[img_pixel_in] == unique_out[img_pixel_out])

    concat_list = list(unique_in[img_pixel_in]) + list(unique_out[img_pixel_out])

    if len(concat_list) == len(set(concat_list)):
        err = """
        The chosen region, for the uploaded image, does not percolate in the fluid flow direction.
        Redo your simulation with an image that makes sense with a proper set of inputs parameters.        
        """
        raise EnvironmentError(err)

    background = np.intersect1d(unique_in, unique_out)

    img[np.isin(labels, background, invert=True)] = no_flow_pixel
    # skimage.io.imsave('/home/rafael/Desktop/am8-2.png', img)

    return img
