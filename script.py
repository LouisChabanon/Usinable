"""
Ce script verifie qu'une surface f est usinable en roulant.
Pour cela on execute la fonction est_usinable(f)
"""

import matplotlib.pyplot as plt
import numpy as np

# Rayon de la fraise
rayon_fraise = 4e-3

# Intervalles de valeurs des variables u et v
intervalle_u = np.arange(-np.pi, np.pi, 0.1)
intervalle_v = np.arange(-np.pi, np.pi, 0.1)

# Fonction f parametrant la surface à étudier


def f(u: float, v: float) -> tuple:
    return v*np.cos(u), v*np.sin(u), v*np.cos(u)*np.sin(u)


"""

Fonctions de derivation

"""


def derivee1_v(f, h=0.001) -> list:
    '''
    Calcule la dérivée partielle de f par rapport à v en tout point
    :param f: Fonction parametrant la surface
    :param h: Petite valeur pour le calcul de dérivée
    :return: Liste de matrice representant les dérivées dans les directions x, y, z
    '''

    # Initialise les matrices pour stocker les dérivées dans les directions x, y, z
    derivees = [np.zeros((len(intervalle_u)-1, len(intervalle_v)-1))
                for _ in range(3)]

    for i in range(3):
        for u in range(len(intervalle_u)-1):
            for v in range(len(intervalle_v)-1):
                derivees[i][u][v] = (
                    f(intervalle_u[u], intervalle_v[v]+h)[i]-f(intervalle_u[u], intervalle_v[v])[i])/h
    return derivees


def derivee1_u(f, h=0.001) -> list:
    '''
    Calcule la dérivée partielle de f par rapport à u en tout point
    :param f: Fonction parametrant la surface
    :param h: Petite valeur pour le calcul de dérivée
    :return: Liste de matrice representant les dérivées dans les directions x, y, z
    '''

    # Initialise les matrices pour stocker les dérivées dans les directions x, y, z
    derivees = np.array(
        [np.zeros((len(intervalle_u)-1, len(intervalle_v)-1)) for _ in range(3)])

    for i in range(3):
        for u in range(len(intervalle_u)-1):
            for v in range(len(intervalle_v)-1):
                derivees[i, u, v] = (
                    f(intervalle_u[u]+h, intervalle_v[v])[i]-f(intervalle_u[u], intervalle_v[v])[i])/h
    return derivees


def derivee2_u(f, h=0.001) -> list:
    '''
    Dérivée partielle seconde discrète par rapport à u, renvoie une liste de 3 matrices
    :param f: Fonction parametrant la surface
    :param h: Petit element de dérivation
    :return: Liste de matrices de dérivées dans les directions x,y,z
    '''

    # Initialise les matrices pour stocker les dérivées secondes de f dans les direction x,y,z
    derivees = np.array(
        [(np.zeros((len(intervalle_u)-1, len(intervalle_v)-1))) for _ in range(3)])
    for i in range(3):
        for u in range(1,len(intervalle_u)-1):
            for v in range(1,len(intervalle_v)-1):
                derivees[i][u][v] = (f(intervalle_u[u]+h, intervalle_v[v])[i]-2*f(
                    intervalle_u[u], intervalle_v[v])[i] + f(intervalle_u[u]-h, intervalle_v[v])[i])/h**2
    return derivees


def derivee2_v(f, h=0.001) -> list:
    '''
    Dérivée partielle seconde discrète par rapport à v, renvoie une liste de 3 matrices
    :param f: Fonction parametrant la surface
    :param h: Petit element de dérivation
    :return: Liste de matrices de dérivées dans les directions x,y,z
    '''

    # Initialise les matrices pour stocker les dérivées secondes de f dans les direction x,y,z
    derivees = [np.zeros((len(intervalle_u)-1, len(intervalle_v)-1))
                for _ in range(3)]
    for i in range(3):
        for u in range(1,len(intervalle_u)-1):
            for v in range(1,len(intervalle_v)-1):
                derivees[i][u][v] = (f(intervalle_u[u], intervalle_v[v]+h)[i]-2*f(
                    intervalle_u[u], intervalle_v[v])[i] + f(intervalle_u[u], intervalle_v[v]-h)[i])/h**2
    return derivees


"""

Fonctions de norme et produit vectoriel

"""


def norme(v: list) -> float:
    '''Renvoie la norme d'un vecteur v'''
    return (v[0]**2 + v[1]**2 + v[2]**2)**0.5


def produit_vectoriel(a: list, b: list) -> list:
    """
    Calcul le produit vectoriel de deux vecteurs a, b

    :param a: Vecteur a
    :param b: Vecteur b
    :return: Produit vectoriel c = a x b
    """
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]
    return c


"""

Fonction de rayon de courbure

"""


def rayon_courbure(f) -> np.ndarray:
    '''
    Calcule le rayon de courbure de f dans la direction u
    :param f: Fonction parametrant la surface
    :return: Matrice du rayon de courbure en chaque point
    '''

    du = derivee1_u(f)
    d2u = derivee2_u(f)
    # la liste des rayons en tout points
    rayon_courbure = np.zeros((len(intervalle_u), len(intervalle_v)))
    for i in range(len(intervalle_u)-1):
        for j in range(len(intervalle_v)-1):
            # On calcule 1/R
            rayon = norme(produit_vectoriel(
                du[:, i, j], d2u[:, i, j])) / norme(du[:, i, j])**3

            # On evite de diviser par zero
            if abs(rayon) < 1e-10:
                # Si 1/R est nul R tend vers l'infini
                rayon_courbure[i][j] = float('inf')
            else:
                rayon_courbure[i][j] = 1/rayon
    return rayon_courbure


def est_usinable(f) -> bool:
    '''
    La fonction verifie que la surface définie par f est bien usinable en roulant
    Pour cela on verifie dans un premier temps que la surface n'est pas réglée (alors la surface n'est pas usinable)
    puis on verifie que la surface à un rayon de courbure superieur à celui de la fraise
    :param f: fonction parametrant la surface
    :Return: True si la surface est usinable
    '''

    d1 = derivee1_v(f)
    l = derivee2_v(f)

    # On verifie que les dérivés secondes par rapport à v sont nulles
    for i in range(0, 3):
        for j in range(len(l[1])):
            for k in range(len(l[1])):
                if l[i][j][k] > 10**(-5):
                    print("La surface est non réglée")
                    return False

    # On cherche maintenant à verifier que le rayon de courbure est superieur à celui de la fraise
    if np.any(rayon_courbure(f)) <= rayon_fraise:
        print("Le rayon de courbure de la courbe est trop faible")
        return False

    print("La surface est usinable !")
    return True


if __name__ == "__main__":
    print(est_usinable(f))
