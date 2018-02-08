from . import IEMD as iemd
import math
import numpy as np

def image_freq_and_ampl(img, eng):
    maxima, minima, maxima_loc, minima_loc = iemd.extrema(img, eng, depth=0)
    amplitude = (np.median(maxima[maxima_loc.astype(bool)]) - np.median(minima[minima_loc.astype(bool)])) / 2
    num_peaks = np.sum(maxima_loc)
    area = img.shape[0] * img.shape[1]
    frequency = math.sqrt(num_peaks / area)
    return frequency, amplitude


def create_masking_signal(img, freq, ampl, freq_modifier, ampl_modifier=1):
    if (freq_modifier < 1) or (freq_modifier > 1.49):
        raise Exception("Frequency modifier must be larger than 1 and smaller than 1.49")
    if freq > 0.5:
        raise Exception("Frequency freq must be lower than 0.5 or else antialising in the masking signal will occur")
    x = np.arange(0, img.shape[0])
    y = np.arange(0, img.shape[1])
    am = ampl * ampl_modifier
    fm = freq * freq_modifier
    xx, yy = np.meshgrid(x, y)
    masking_signal = am * np.sin(2*np.pi*fm*(yy - 0.5)) + am * np.sin(2*np.pi*fm*(xx - 0.5))
    return masking_signal


def mode_mixing_separation(img, freq, ampl, eng, rbase=15, num_siftings=10, freq_modifier=1.2, ampl_modifier=1, debug=False):
    if debug:
        print("freq_modifier:", freq_modifier, "ampl_modifier:", ampl_modifier)
    masking_signal = create_masking_signal(img, freq, ampl, freq_modifier)
    mask_pos = img + masking_signal
    mask_neg = img - masking_signal
    imf_pos, residue_pos = iemd.IEMD(mask_pos, eng, rbase=rbase, num_siftings=num_siftings, max_imfs=1, debug=debug)
    imf_neg, residue_neg = iemd.IEMD(mask_neg, eng, rbase=rbase, num_siftings=num_siftings, max_imfs=1, debug=debug)
    residue = (residue_pos + residue_neg) / 2
    first_imf = (imf_pos + imf_neg) / 2
    return first_imf, residue

def mode_mixing_separation_noise_removal(img, eng, rbase=15, num_siftings=10, freq_modifier=1.2, ampl_modifier=1, debug=False):
    freq, ampl = image_freq_and_ampl(img, eng)
    return mode_mixing_separation(img, freq, ampl, eng, rbase=rbase, num_siftings=num_siftings, freq_modifier=freq_modifier, debug=debug)
