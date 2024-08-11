import pandas as pd
import cv2
import numpy as np
from skimage.feature import hog, graycomatrix, graycoprops
from scipy.stats import entropy
from skimage.measure import label, regionprops
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from skimage.feature import greycomatrix, greycoprops
from scipy.stats import entropy, skew, kurtosis
from skimage.measure import label, regionprops, moments_hu


def extract_features_from_image(image_path):
    image = cv2.imread(image_path)
    if image is None:
        print(f"Error: Unable to read image {image_path}")
        return None

    original_height, original_width = image.shape[:2]

    # Resize the image
    image_resized = cv2.resize(image, (100, 100))

    # Convert to grayscale for some features
    image_gray = cv2.cvtColor(image_resized, cv2.COLOR_RGB2GRAY)
    image_resized_gray = cv2.resize(image_gray, (100, 100))


    # HOG features
    hog_features = hog(image_gray, orientations=8, pixels_per_cell=(8, 8),
                       cells_per_block=(2, 2), visualize=False)
    hog_mean = np.mean(hog_features)
    hog_std = np.std(hog_features)


    # Texture features from grayscale
    glcm = graycomatrix(image_resized_gray, [5], [0], symmetric=True, normed=True)
    contrast = graycoprops(glcm, 'contrast')[0, 0]
    dissimilarity = graycoprops(glcm, 'dissimilarity')[0, 0]
    homogeneity = graycoprops(glcm, 'homogeneity')[0, 0]
    ASM = graycoprops(glcm, 'ASM')[0, 0]
    energy = graycoprops(glcm, 'energy')[0, 0]
    correlation = graycoprops(glcm, 'correlation')[0, 0]
    entropy_value = entropy(image_resized_gray.ravel())
    mean_intensity = np.mean(image_resized_gray)
    std_intensity = np.std(image_resized_gray)
    skewness_intensity = skew(image_resized_gray.ravel())
    kurtosis_intensity = kurtosis(image_resized_gray.ravel())

    # Shape features
    edges = cv2.Canny(image_gray, 100, 200)
    num_centromeres = len(regionprops(label(edges)))
    chrom_length = original_height - 2
    aspect_ratio = original_height / original_width

    # Hu Moments
    moments = moments_hu(image_gray)
    hu_moments = [m for m in moments]

    # Edge-based features
    # edge_density = np.sum() / (original_height * original_width)
    _, binary_image = cv2.threshold(image_resized_gray, 20, 255, cv2.THRESH_BINARY)

    # Contour-based features
    contours, _ = cv2.findContours(binary_image, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    total_edge_pixels = np.sum([cv2.contourArea(contour) for contour in contours])
    edge_density = total_edge_pixels / (original_height * original_width)
    if contours:
        largest_contour = max(contours, key=cv2.contourArea)
        perimeter_length = cv2.arcLength(largest_contour, True)
        solidity = cv2.contourArea(largest_contour) / cv2.contourArea(cv2.convexHull(largest_contour))
    else:
        perimeter_length = 1
        solidity = 0


    # Color features
    mean_color = np.mean(image_resized, axis=(0, 1))
    std_color = np.std(image_resized, axis=(0, 1))

    features = {
        'hog_mean': hog_mean,
        'hog_std': hog_std,
        'dissimilarity': dissimilarity,
        'homogeneity': homogeneity,
        'ASM': ASM,
        'energy': energy,
        'correlation': correlation,
        'entropy': entropy_value,
        'mean_intensity': mean_intensity,
        'std_intensity': std_intensity,
        'skewness_intensity': skewness_intensity,
        'kurtosis_intensity': kurtosis_intensity,
        'aspect_ratio': aspect_ratio,
        'edge_density': edge_density,
        'solidity': solidity,
        'mean_color_R': mean_color[0],
        'std_color_R': std_color[0],
    }

    return features