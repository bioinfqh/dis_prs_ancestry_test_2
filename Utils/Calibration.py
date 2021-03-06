import numpy as np
import matplotlib.pyplot as plt
import calibration as cal
from sklearn.isotonic import IsotonicRegression
from sklearn import preprocessing
from sklearn.calibration import  calibration_curve
from sklearn.metrics import confusion_matrix, plot_confusion_matrix,brier_score_loss,\
mean_squared_error,log_loss,precision_score, f1_score, recall_score, roc_auc_score

"""
This code uses library uncertianity-calibration from 
@inproceedings{kumar2019calibration,
  author = {Ananya Kumar and Percy Liang and Tengyu Ma},
  booktitle = {Advances in Neural Information Processing Systems (NeurIPS)},
  title = {Verified Uncertainty Calibration},
  year = {2019},
}
"""

def check_prob_sum(probs):
    """
    Check that the vector probs sums to 1
    """
    rtol = 0.1
    try:
        np.testing.assert_allclose(sum(probs),1, rtol = rtol)
    except AssertionError:
        print(f'AssertionError: probs do not sum to within {rtol} of 1')



def normalize_prob(prob,n_classes):
    """ Normalize the probabilities
    """
    if n_classes == 2:
        prob[:, 0] = 1. - prob[:, 1]
    else:
        prob /= np.sum(prob, axis=1)[:, np.newaxis]

    # XXX : for some reason all probas can be 0
    prob[np.isnan(prob)] = 1. / n_classes

    # Deal with cases where the predicted probability minimally exceeds 1.0
    prob[(1.0 < prob) & (prob <= 1.0 + 1e-5)] = 1.0
    
    return prob


def plot_reliability_curve(pred_prob,y_cal,pop_order,method = 'Uncalibrated', bins = 10):
    
    fig, (ax1) = plt.subplots(nrows = 1, ncols=1, figsize = (8,6))
                                              
    ax1.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")
    ax1.set_title(method + ' Reliability Plot')
        
    class_metrics = {}
    n_classes = len(pop_order)
    for i in range(n_classes):
        y_label = np.where(y_cal==i,1,0)
        est_prob = pred_prob[:,:,i]
        true_prob, pred_proba = calibration_curve(y_label.flatten(),est_prob.flatten(), n_bins = 15)
        ax1.plot(pred_proba, true_prob, label = pop_order[i])
        
        ax1.set_xlabel('Est. Prob/ mean predicted value')
        ax1.set_ylabel('True Prob/fraction_of_positives')
        ax1.legend(loc="lower right")
        
        y_test = y_label.flatten()
        y_pred =np.where(est_prob.flatten()>1/n_classes ,1,0)
        
        class_metrics[pop_order[i]]={}
        class_metrics[pop_order[i]]['Precision'] = float("%1.3f" %precision_score(y_test, y_pred))
        class_metrics[pop_order[i]]['Recall'] = float("%1.3f" %recall_score(y_test, y_pred))
        class_metrics[pop_order[i]]['F1'] = float("%1.3f" %f1_score(y_test, y_pred))
        class_metrics[pop_order[i]]['AUC'] = float("%1.3f" %roc_auc_score(y_test, y_pred))

    
    fig.tight_layout()
    plt.show()
    return class_metrics


def calibrator_module(zs,ys,n_classes,method):
        
    if method == 'Platt':
        calibrator = cal.PlattBinnerMarginalCalibrator(n_classes, num_bins=200)
        calibrator.train_calibration(zs, ys)
      
        calibrated_zs = calibrator.calibrate(zs)
        calibration_error = cal.get_calibration_error(calibrated_zs, ys)
        print("Scaling-binning L2 calibration error with %s is %.2f%%" % (method, (100 * calibration_error)))
        
    elif method =='Isotonic':

        # binarize class labels
        lb = preprocessing.LabelBinarizer().fit(ys)
        y_cal_ohe = lb.transform(ys)
        del lb
        
        calibrator = []
        # calibrated_prob = np.zeros((ys.shape[0],n_classes))
        
        for i in range(n_classes):    
            calibrator.append(IsotonicRegression(out_of_bounds = 'clip').fit(zs[:,i], y_cal_ohe[:,i]))
            # calibrated_prob[:,i] = calibrator[i].transform(zs[:,i])

        # # Normalize the probabilities
        # calibrated_prob = normalize_prob(calibrated_prob,n_classes)
        # calibration_error = cal.get_calibration_error(calibrated_prob, ys)
        # print("Scaling-binning L2 calibration error with %s is %.2f%%" % (method, (100 * calibration_error)))
        
    return calibrator


def comparison(uncalibrated_zs,calibrated_zs,ys,pop_type,pop_order):
    
    """
    Compare probabilities of a population type before and after calibration
    """
    # pop_type (int between 0 and len(pop_order) for which population the comparison needs to be
    
    fig, (ax1,ax2) = plt.subplots(nrows = 2, ncols=1, figsize = (10,10))
    
    y_label = np.where(ys==pop_type,1,0)
    uncal_est_prob = uncalibrated_zs[:,:,pop_type]
    cal_est_prob = calibrated_zs[:,:,pop_type]
    uncal_true_prob, uncal_pred_proba = calibration_curve(y_label.flatten(),uncal_est_prob.flatten(), n_bins = 10)
    ax1.plot(uncal_pred_proba, uncal_true_prob, label = pop_order[pop_type] + ' uncalib.')
    
    cal_true_prob, cal_pred_proba = calibration_curve(y_label.flatten(),cal_est_prob.flatten(), n_bins = 10)
    ax1.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")
    ax1.plot(cal_pred_proba, cal_true_prob, label = pop_order[pop_type] + ' calib.')
    
    ax1.set_xlabel('Est. Prob/ mean predicted value')
    ax1.set_ylabel('True Prob/fraction_of_positives')
    ax1.legend(loc="lower right")
    
    ax2.hist(uncal_est_prob.flatten(), range=(0, 1), bins=20, histtype="step", label = 'uncalib.', log=True, lw=3)
    ax2.hist(cal_est_prob.flatten(), range=(0, 1), bins=20, histtype="step", label = 'calib.', log=True, lw=3)
    ax2.legend(loc="lower right")
        
    calibration_error = cal.get_ece(uncal_est_prob.flatten(), y_label.flatten())
    print("Scaling-binning L2 calibration error for uncalib. probs is %.2f%%" % (100 * calibration_error))
    
    calibration_error = cal.get_ece(cal_est_prob.flatten(), y_label.flatten())
    print("Scaling-binning L2 calibration error for calib. probs is %.2f%%" % (100 * calibration_error))
        
    fig.tight_layout()
    plt.show()