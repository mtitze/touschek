import matplotlib.pyplot as plt
import numpy as np

def make_survey_labels(x, y, keywords, thetas):
    no_drift = np.where(np.logical_and(keywords != 'drift', keywords != 'marker')) # do not label drifts
    x = x[no_drift]
    y = y[no_drift]
    keywords = keywords[no_drift]
    thetas = thetas[no_drift]
    
    n_elements = len(keywords)
    x_label, y_label = [], []
    labels = []
    thetas_new = []
    
    current_label = keywords[0]
    current_indices = [0]
    
    for k in range(1, n_elements):
        keyword = keywords[k]
        if keyword == current_label:
            current_indices.append(k)
            continue
        else:
            # a group of points and labels has been determined, given by the indices in 'current_indices'.
            nn = int(len(current_indices)/2)
            if nn > 0:
                index = current_indices[nn]
            else:
                index = k - 1
            x_label.append(x[index])
            y_label.append(y[index])
            labels.append(keywords[index])
            thetas_new.append(thetas[index])
            current_indices = [k]
            current_label = keyword
            
    return x_label, y_label, labels, thetas_new


def plot_survey(madx, kmin=None, kmax=None, figsize=(12, 12), s=3, aspect=True):
    '''
    Plot survey of machine elements utilizing MAD-X survey command.
    
    INPUT
    =====
    aspect: If set to True, create plot with same aspect ratio in both directions.
    '''
    surv = madx.survey()

    x = surv.x[kmin:kmax]
    z = surv.z[kmin:kmax]
    labels = surv.keyword[kmin:kmax]
    thetas = surv.theta[kmin:kmax]

    xl, zl, labelsl, thetasl = make_survey_labels(x, z, labels, thetas)

    plt.figure(figsize=figsize)
    plt.scatter(x, z, s=s)

    for i, txt in enumerate(labelsl):
        angle = -thetasl[i]*180/np.pi
        if angle > 90 and angle < 270:
            text = txt + ' '*(2*len(txt) + 1)
            angle = angle + 180
        else:
            text = ' '*(2*len(txt) + 1) + txt
        plt.text(xl[i], zl[i], text, rotation=angle, ha='center', va='center')
    plt.xlabel(r'$x$ [m]')
    plt.ylabel(r'$z$ [m]')
    ax = plt.gca()
    ax.set_aspect(1)
    return plt

def plot_touschek_losses(optics, touschek_results, xlim=[], with_beta=False,
                         figsize=(16, 4)):

    plt.figure(figsize=figsize)
    plt.title(f'Touschek lifetime: {touschek_results["lifetime"]/60/60:.3f} [h]', loc='right')
    pos = optics.function.position.values
    values = touschek_results['touschek_const']*touschek_results['touschek_ring']    
    if len(xlim) > 0:
        lim_indices = np.logical_and(pos >= xlim[0], pos <= xlim[-1])
        pos = pos[lim_indices]
        values = values[lim_indices]

    plt.plot(pos, values, label='losses')
    plt.ylabel(r'$\frac{r_p^2 c N_p F(\tau_m, B_1, B_2)}{8 \pi \gamma^2 \tau_m \sigma_s \sqrt{\sigma_x^2 \sigma_y^2 - \delta^4 D_x^2 D_y^2}}$ [1/s]',
            fontsize=14)


    exclude_in_ticks = ['marker']
    keywords = optics.madx.table.twiss.keyword
    s = optics.madx.table.twiss.s
    if len(xlim) > 0:
        lim_indices2 = np.logical_and(s >= xlim[0], s <= xlim[-1])
        s = s[lim_indices2]
        keywords = keywords[lim_indices2]
    indices_oi = [k for k in range(len(keywords)) if keywords[k] not in exclude_in_ticks]
    plt.xticks(s[indices_oi], keywords[indices_oi], rotation=90)
    plt.twiny()
    plt.scatter([pos[0], pos[-1]], [values[-1], values[-1]], alpha=0)
    plt.xlabel(r'$s$ [m]', fontsize=14)

    if with_beta:
        plt.twinx()
        bxbyi = 1/(optics.function.betax*optics.function.betay)
        if len(xlim) > 0:
            bxbyi = bxbyi[lim_indices]
        plt.plot(pos, bxbyi, color='black', linestyle='--', label=r'$1/(\beta_x \beta_y)$')
        plt.ylabel(r'$1/(\beta_x \beta_y)$', fontsize=14)
        plt.legend()

    plt.show()
