from numpy.core.fromnumeric import shape
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

def plot_algorithm_performance(c_before_after_df, dmc_before_after_df, func, dim, ax): 
    #f, ax = plot.subplots(figsize=fig_size)
   
    t=c_before_after_df.query("Func == '" + func + "' and Dim == " + str(dim))
    x=pd.DataFrame([t.iloc[0].index[2:],t.iloc[0].values[2:]]).T

    yt=dmc_before_after_df.query("Func == '" + func + "' and Dim == " + str(dim))
    y=pd.DataFrame([yt.iloc[0].index[2:],yt.iloc[0].values[2:]]).T
    
    x.plot(ax = ax)
    y.plot(ax = ax)

    ax.set_title(func + " with D=" + str(dim))
    ax.legend(['C-GRASP', 'DMC-GRASP'])
    #plot.show()

def plot_alg_perf_by_seed(seed):
	c_before_after_df = pd.read_csv("../BEFORE_AFTER/C/before_after_"+ str(seed) +".csv", sep=";")
	dmc_before_after_df = pd.read_csv("../BEFORE_AFTER/DMC/before_after_"+ str(seed) +".csv", sep=";")

	list_of_instances = [('BOOTH', 2), ('BRANIN', 2), ('EASOM', 2), ('GOLDSTEINPRICE', 2), ('MATYAS', 2), ('HUMP', 2), ('ROSENBROCK', 2), ('SCHWEFEL', 2), ('SHUBERT', 2), ('ZAKHAROV', 2), ('SPHERE', 3), ('HARTMANN', 3), ('COLVILLE', 4), ('PERM', 4), ('PERM0', 4), ('POWERSUM', 4), ('SHEKEL', 5), ('SHEKEL', 7), ('SHEKEL', 10), ('ROSENBROCK', 5), ('ZAKHAROV', 5), ('HARTMANN', 6), ('SCHWEFEL', 6), ('TRID', 6), ('GRIEWANK', 10), ('RASTRIGIN', 10), ('ROSENBROCK', 10), ('SUMSQUARES', 10), ('TRID', 10), ('ZAKHAROV', 10), ('GRIEWANK', 20), ('RASTRIGIN', 20), ('ROSENBROCK', 20), ('SUMSQUARES', 20), ('ZAKHAROV', 20), ('POWELL', 24), ('DIXONPRICE', 25), ('ACKLEY', 30), ('LEVY', 30), ('SPHERE', 30)]
	for instance1, instance2 in zip(list_of_instances[0::2], list_of_instances[1::2]):
		fig, (ax1, ax2) = plt.subplots(1, 2, figsize=fig_size)
		plot_algorithm_performance(c_before_after_df, dmc_before_after_df, instance1[0], instance1[1], ax1)
		plot_algorithm_performance(c_before_after_df, dmc_before_after_df, instance2[0], instance2[1], ax2)
		st.pyplot(plt)


st.title('Análise exploratória dos resultados')
fig_size = (12, 3)

window = st.slider("Sementes", min_value = 1)
seed = 270000 + window
plot_alg_perf_by_seed(seed)


