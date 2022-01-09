"""
Computes channel capacity using the Blahut-Arimoto algorithm.
"""
import sys
import numpy as np
import math

def get_CC(_response):
    r
    """
    Computes channel capacity using the Blahut-Arimoto algorithm.

    Parameters
    ----------
    _response: numpy.array
        2D numpy array, with number of rows as the number of input levels for the channel, and the number of columns is equal
        the number of output levels. This array is the transition matrix for the information channel.

    Returns
    -------
    C: float
        Channel capacity of the information channel in bits.

    err: float
        Error value at which the iterative Blahut-Arimoto algorithm terminates.

    prob: numpy.array
        1D numpy.array with the optimal input distribution that achieves the channel capacity for the information channel.
        The length of this array is equal to the number of rows in _response.
    """

    response = _response

    rows, cols = response.shape[0], response.shape[1]

    # Assume uniform distribution as initial condition
    probs = np.ones(shape=(1,rows))/float(rows)

    # Tolerance for convergence
    errtol = 1e-3

    # Initial value for error to start the iterative loop
    err = 1
    steps = 0

    while err>errtol:
        c = np.zeros(shape=(rows,))

        v1 = np.matmul(probs,response)[0,:]

        for j in range(0,rows):
            for k in range(0,cols):
                if response[j,k]>0.0:
                    c[j] += response[j,k]*math.log(response[j,k]/v1[k])

            c[j] = math.exp(c[j])

        mean_c = np.dot(probs,c)[0]

        I_L = math.log(mean_c)
        I_U = math.log(np.max(c))

        err = abs(I_U - I_L)

        err *= 1.0/math.log(2)

        if err>errtol:
            probs[0,:] = np.multiply(probs[0,:],c)/mean_c
        else:
            C = I_L

        steps += 1

    C *= 1.0/math.log(2)

    return C, err, probs[0]
