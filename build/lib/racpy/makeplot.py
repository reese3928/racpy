
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def makeplot(res, main="RNA age vs chronological age",
             xlab="chronological age", ylab="RNA Age"):
    """Make plot to visualize RNA age.

    This function makes plots to visualize the relationship between
    chronological age and RNA age.

    :param res: a pandas DataFrame returned by predict_age function.
    If the chronological age is not provided when using predict_age
    function, visulization cannot be made.

    :param main: a string which specifies the title of the plot
    :param xlab: a string which specifies label of x-axis
    :param ylab: label of y-axis
    :return: a scatter plot
    """
    assert isinstance(res, pd.DataFrame), \
        "res should be a pandas DataFrame."
    assert "ChronAge" in res.columns, \
        "Chronological age is not found in res DataFrame."
    if res.dropna().shape[0] < 30:
        print("Less than 30 samples. The linear regression on RNA age vs "
              "chronological age may not be reliable.")

    myplot = sns.lmplot(x='ChronAge', y='RNAAge', data=res, fit_reg=True)
    myplot.set_axis_labels(xlab, ylab)
    plt.title(main)
    return myplot
