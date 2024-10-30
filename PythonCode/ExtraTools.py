import numpy as np
import matplotlib.pyplot as plt


def metricsComputation(pricesLewis, pricesModel):
    normInf = np.max(np.abs(pricesLewis - pricesModel))
    RMSE = np.sqrt(np.mean((pricesLewis - pricesModel) ** 2))
    MAPE = np.mean(np.abs((pricesLewis - pricesModel) / pricesLewis)) * 100

    return normInf, RMSE, MAPE


def plotPrices(moneyness, pricesLewis, prices_ED, conf_int_ED, prices_FM, conf_int_FM, plot_title):
    plt.figure(figsize=(12, 6))
    plt.plot(moneyness, pricesLewis, label='Lewis FFT', color='blue')
    plt.plot(moneyness, prices_ED, label='Exact Decomposition', color='green')
    plt.plot(moneyness, conf_int_ED[0], label='Lower bound Exact Decomposition', color='green', alpha=0.4)
    plt.plot(moneyness, conf_int_ED[1],  label='Upper bound Exact Decomposition', color='green', alpha=0.4)
    plt.plot(moneyness, prices_FM, label='Fast General Monte Carlo', color='red')
    plt.plot(moneyness, conf_int_FM[0], label='Lower bound FGMC', color='red', alpha=0.4)
    plt.plot(moneyness, conf_int_FM[1], label='Upper bound FGMC', color = 'red', alpha=0.4)

    plt.xlabel('Moneyness')
    plt.ylabel('Prices')
    plt.title( plot_title)
    plt.legend()
    plt.grid(True)
    plt.show()

