import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

def plot_RDI(data):
    indices = ['PI', 'DI', 'PDI']
    timepoints = ['0H', '1H', '2H', '3H', '4H']

    data['TimePoint'] = pd.Categorical(data['TimePoint'], categories=timepoints, ordered=True)

    for index in indices:
        for cytokine in data['Cytokines'].unique():
            plt.figure(figsize=(8, 6))
            data_cyt = data[data['Cytokines'] == cytokine]
            ax = sns.boxplot(data=data_cyt, x='TimePoint', y=index, hue='Type', showfliers=False,
                             palette=sns.color_palette('pastel'))
            timepoint_positions = {t: i for i, t in enumerate(timepoints)}  

            for timepoint in timepoints:
                data_tp = data_cyt[data_cyt['TimePoint'] == timepoint]

                types = data_tp['Type'].unique()
                for i in range(len(types)):
                    if i != (len(types) - 1):
                        x = data_tp[data_tp['Type'] == types[i]][index].values
                        y = data_tp[data_tp['Type'] == types[i+1]][index].values

                    if i == (len(types) - 1):
                        x = data_tp[data_tp['Type'] == types[i]][index].values
                        y = data_tp[data_tp['Type'] == types[0]][index].values


                    if len(x) > 0 and len(y) > 0:  # Alleen testen als beide groepen data hebben
                        stat, p_value = mannwhitneyu(x, y, alternative='two-sided')
                        print(p_value)
                        # 🔥 Betekenisvolheid aangeven met sterretjes
                        significance = 'n.s.'
                        if p_value < 0.05: significance = '*'
                        if p_value < 0.01: significance = '**'
                        if p_value < 0.001: significance = '***'

                        x_pos = timepoint_positions[timepoint]
                        y_pos = 2

                        if i != (len(types)-1):
                            fraction = 0.8/len(types)
                            x_pos_left = x_pos - 0.4 + (i*fraction) + (fraction/2)
                            x_pos_right = x_pos - 0.4 + ((i+1) * fraction) + (fraction / 2)

                            plt.text((x_pos_left+x_pos_right)/2, y_pos, significance, ha='center', fontsize=12, fontweight='bold', color='red')
                            plt.plot([x_pos_left, x_pos_right], [y_pos - 0.05, y_pos - 0.05], color='black', lw=1.5)

                        if i == (len(types) - 1):
                            fraction = 0.8 / len(types)
                            x_pos_left = x_pos - 0.4 + (i * fraction) + (fraction / 2)
                            x_pos_right = x_pos - 0.4 + (fraction / 2)

                            plt.text((x_pos_left + x_pos_right) / 2, y_pos + 0.15, significance, ha='center', fontsize=12,
                                     fontweight='bold', color='red')
                            plt.plot([x_pos_left, x_pos_right], [y_pos + 0.1, y_pos + 0.1], color='black', lw=1.5)

                    if len(types) == 2:
                        break

            plt.title(f'{index} for cytokine: {cytokine}')
            plt.ylim(0, 2.3)
            plt.legend(title='Type', bbox_to_anchor = (1.1, 1.1), loc=0)
            plt.show()

    return
