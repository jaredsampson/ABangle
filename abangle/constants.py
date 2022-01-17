from typing import List

coresets = {
    'H': {
        35, 12, 38, 36, 83, 19, 94, 
        37, 11, 47, 39, 93, 46, 45, 
        68, 69, 71, 70, 17, 72, 92, 
        84, 91, 90, 20, 21, 85, 25, 
        24, 86, 89, 88, 87, 22, 23
    },
    'L': {
        44, 19, 69, 14, 75, 82, 15, 
        21, 47, 20, 48, 49, 22, 81, 
        79, 80, 23, 36, 35, 37, 74,
        88, 38, 18, 87, 17, 86, 85,
        46, 70, 45, 16, 71, 72, 73
    }
}

principle_components = {
    'H': [
        [0.9525187, -0.1371821, 0.2718256],
        [-0.117058, 0.659152, 0.7428432],
        [-2.691829, -3.847092, 1.196887]
    ],
    'L': [
        [-0.6193343, 0.639472, 0.4555223],
        [0.5267385, 0.7686645, -0.362907],
        [-3.702842, -0.6288583, -5.314558]
    ]
}