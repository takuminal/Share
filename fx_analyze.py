import matplotlib.pyplot as plt
import yfinance as yf
import numpy as np
import pandas as pd

ticker = "JPYUSD=X"

tick = yf.Ticker(ticker)
# get historical market data
hist = tick.history(period = "1d", interval="1m")

change = hist["Close"] - hist["Open"]
t = change.index
val = change.values

l = len(val)

split_num= 100
t_hour = np.array_split(t.hour,split_num)
t_min = np.array_split(t.minute,split_num)

t_head = [f"{str(t_hour[i][0])}:{str(t_min[i][0])}" for i in range(split_num)]
val_split = np.array_split(val, split_num)

sum_arr = []
mean_arr = []

val_sum = 0

for elem in val_split:
    val_sum = elem.sum() + val_sum
    val_mean =elem.mean()
    sum_arr.append(val_sum)
    mean_arr.append(val_mean)
    
plt.figure(figsize=(10,5))
plt.subplot(2, 1, 1) 
plt.boxplot(val_split,positions=[i for i in range(split_num)])
plt.plot(t_head,mean_arr)
plt.xticks(range(0, split_num),t_head,rotation=90)
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(t_head,sum_arr)
plt.xticks(range(0, split_num),t_head,rotation=90)
plt.grid(True)



plt.show()