# This file is part of Frhodo. Copyright © 2020, UChicago Argonne, LLC
# and licensed under BSD-3-Clause. See License.txt in the top-level 
# directory for license and copyright information.

import colors

import matplotlib as mpl
import numpy as np
from scipy import stats

from plot.base_plot import Base_Plot


colormap = colors.colormap(reorder_from=1, num_shift=4)

class Plot(Base_Plot):
    def __init__(self, parent, widget, mpl_layout):
        super().__init__(parent, widget, mpl_layout)
        self.canvas.mpl_connect("motion_notify_event", self.hover)
        
    def create_canvas(self):
        self.ax = []
        
        ## Set left (QQ) plot ##
        self.ax.append(self.fig.add_subplot(1,2,1))
        self.ax[0].item = {}
        self.ax[0].item['qq_data'] = []
        self.ax[0].item['ref_line'] = self.ax[0].add_line(mpl.lines.Line2D([],[], c=colormap[0], zorder=2))
        
        self.ax[0].text(.5,1.03,'QQ-Plot of Residuals', fontsize='large',
            horizontalalignment='center', verticalalignment='top', transform=self.ax[0].transAxes)
        
        self.fig.subplots_adjust(left=0.06, bottom=0.065, right=0.98,
                        top=0.965, hspace=0, wspace=0.12)
        
        ## Set right (density) plot ##
        self.ax.append(self.fig.add_subplot(1,2,2))
        self.ax[1].item = {}
        
        self.ax[1].item['density'] = [self.ax[1].add_line(mpl.lines.Line2D([],[],
                                                          ls='-',  c=colormap[0], zorder=3))]
        
        # add exp scatter/line plots
        self.add_exp_plots()
        self.ax[1].item['shade'] = [self.ax[1].fill_between([0, 0], 0, 0)]
        
        self.ax[1].item['annot'] = self.ax[1].annotate("", xy=(0,0), xytext=(-100,20), 
                        textcoords="offset points", bbox=dict(boxstyle="round", fc="w"), 
                        arrowprops=dict(arrowstyle="->"), zorder=10)
        self.ax[1].item['annot'].set_visible(False)
        
        self.ax[1].item['outlier_l'] = []
        for i in range(2):
            self.ax[1].item['outlier_l'].append(self.ax[1].axvline(x=np.nan, ls='--', c='k'))# c='#BF0000'))
        
        self.ax[1].text(.5,1.03,'Density Plot of Residuals', fontsize='large',
            horizontalalignment='center', verticalalignment='top', transform=self.ax[1].transAxes) 

        # Create canvas from Base
        super().create_canvas()
     
    def add_exp_plots(self):
        i = len(self.ax[1].item['density'])
        
        line = self.ax[1].add_line(mpl.lines.Line2D([],[], color=colormap[i+1], alpha=0.5, zorder=2))
        line.set_pickradius(line.get_pickradius()*2)
        self.ax[1].item['density'].append(line)
        
        scatter = self.ax[0].scatter([],[], color=colormap[i+1], facecolors=colormap[i+1], 
            s=16, linewidth=0.5, alpha = 0.85)
        scatter.set_pickradius(scatter.get_pickradius()*2)
        self.ax[0].item['qq_data'].append(scatter)
     
    def hover(self, event):
        def update_annot(line, ind):
            def closest_xy(point, points):
                dist_sqr = np.sum((points - point[:, np.newaxis])**2, axis=1)
                return points[:,np.argmin(dist_sqr)]
                
            annot = self.ax[1].item['annot']
            x,y = line.get_data()
            xy_mouse = self.ax[1].transData.inverted().transform([event.x, event.y])
            annot.xy = closest_xy(xy_mouse, np.array(line.get_data())[:, ind['ind']])   # nearest point to mouse
            
            text = '{:s}\nExp # {:d}'.format(line.shock_info['series_name'], line.shock_info['num'])
            annot.set_text(text)
            extents = annot.get_bbox_patch().get_extents()
            if np.mean(self.ax[1].get_xlim()) < annot.xy[0]: # if on left side of plot
                annot.set_x(-extents.width*0.8+20)
            else:
                annot.set_x(0)
            annot.get_bbox_patch().set_alpha(0.85)
            annot.set_visible(True)
    
        # if the event happened within axis and no toolbar buttons active do nothing
        if event.inaxes != self.ax[1] or self.toolbar._active is not None: return
        
        draw = False    # tells to draw or not based on annotation visibility change
        default_pick_radius = self.ax[1].item['density'][1].get_pickradius()
        contains_line = []
        for line in self.ax[1].item['density'][1:]:
            contains, ind = line.contains(event)
            if contains and hasattr(line, 'shock_info'):
                contains_line.append(line)     
        
        if len(contains_line) > 0:      # reduce pick radius until only 1 line contains event
            for n in range(1, default_pick_radius)[::-1]:
                if len(contains_line) == 1: # if only 1 item in list break
                    break
                    
                for i in range(len(contains_line))[::-1]:
                    if len(contains_line) == 1: # if only 1 item in list break
                        break
                    
                    contains_line[i].set_pickradius(n)
                    contains, ind = contains_line[i].contains(event)
                    
                    if not contains:
                        del contains_line[i]
            
            # update annotation based on leftover
            contains, ind = contains_line[0].contains(event)
            update_annot(contains_line[0], ind)
            draw = True
            
            for line in contains_line:  # resset pick radius
                line.set_pickradius(default_pick_radius)
                
        elif self.ax[1].item['annot'].get_visible():   # if not over a line, hide annotation
            draw = True
            self.ax[1].item['annot'].set_visible(False)
        
        if draw:
            self.canvas.draw_idle()
     
    def update(self, data, update_lim=True):
        def shape_data(x, y): return np.transpose(np.vstack((x, y)))
        
        shocks2run = data['shocks2run']
        resid = data['resid']
        weights = data['weights']
        resid_outlier = data['resid_outlier']
        num_shocks = len(shocks2run)
        
        # operations needed for both QQ and Density Plot
        allResid = np.concatenate(resid, axis=0)
        # weights = np.concatenate(weights, axis=0)
        # mu = np.average(allResid, weights=weights)
        # allResid -= mu
        # allResid = allResid[np.abs(allResid) < resid_outlier]
        # allResid += mu
        
        fitres = data['fit_result']
        x_grid = np.linspace(fitres[1]-4*fitres[2], fitres[1]+4*fitres[2], 300)
        xlim_density =[x_grid[0], x_grid[-1]]
        
        # add exp line/data if not enough
        for i in range(num_shocks):
            if len(self.ax[1].item['density'])-2 < i:   # add line if fewer than experiments
                self.add_exp_plots()
        
        self.clear_plot(ignore='text', draw=False)
        
        # from timeit import default_timer as timer
        # start = timer()
        
        # for shock in [0]:
            # self.ax[0].add_line(mpl.lines.Line2D([],[], marker='$'+u'\u2195'+'$', 
                # markersize=12, markerfacecolor='#BF0000', markeredgecolor='None', linestyle='None', zorder=4))
        
        # Update left plot
        xrange = np.array([])
        for i in range(num_shocks):
            QQ = data['QQ'][i]
            self.ax[0].item['qq_data'][i].set_offsets(QQ)
            
            xrange = np.append(xrange, [QQ[:,0].min(), QQ[:,0].max()])

            xrange = np.reshape(xrange, (-1,2))
            xrange = [np.min(xrange[:,0]), np.max(xrange[:,1])]
            
        self.ax[0].item['ref_line'].set_xdata(xrange)
        self.ax[0].item['ref_line'].set_ydata(xrange)
        
        # Update right plot
        # clear shades
        for shade in self.ax[1].item['shade'][::-1]:
            shade.remove()
            
        self.ax[1].item['shade'] = []
        
        # kernel density estimates
        gennorm_fit = stats.gennorm.pdf(x_grid, *fitres)
        self.ax[1].item['density'][0].set_xdata(x_grid)
        self.ax[1].item['density'][0].set_ydata(gennorm_fit)
        
        self.ax[1].item['outlier_l'][0].set_xdata(-resid_outlier*np.ones((2,1)))
        self.ax[1].item['outlier_l'][1].set_xdata(resid_outlier*np.ones((2,1)))
        
        for i in range(num_shocks):
            x_grid = data['KDE'][i][:,0]
            density = data['KDE'][i][:,1]
            self.ax[1].item['density'][i+1].set_xdata(x_grid)
            self.ax[1].item['density'][i+1].set_ydata(density)
            self.ax[1].item['density'][i+1].shock_info = shocks2run[i]
                   
            zorder = self.ax[1].item['density'][i+1].zorder
            color = self.ax[1].item['density'][i+1]._color
            shade = self.ax[1].fill_between(x_grid, 0, density, alpha=0.01, zorder=zorder, color=color)
            self.ax[1].item['shade'].append(shade)
        
        if update_lim:
            self.update_xylim(self.ax[0])
            self.update_xylim(self.ax[1], xlim=xlim_density)
    
        # print('{:0.1f} us'.format((timer() - start)*1E3))
            
