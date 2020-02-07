import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from itertools import product
import numpy as np
from scipy.interpolate import splprep,splev

def prepString(shape_str):
    return [int(i) for i in shape_str.split()][2:]

def VisualiseSingleShape(shape,ax=None,corner=(0,1),extent=1,add_direction=False):
    cols=['darkgreen','royalblue','firebrick','goldenrod','mediumorchid']*10
    hatchs=['//','.','\\\\','O']*20
    ar_offsets={0:(0,-.25,0,.25),1:(-.25,0,.25,0),2:(0,.25,0,-.25),3:(.25,0,-.25,0)}
    rescale=False
    if ax is None:
        rescale=True
        _, ax = plt.subplots()

    dx,dy=shape[:2]
    max_d=max(dx,dy)
    extentS=extent/max_d
    for i,j in product(range(dx),range(dy)):
        if(shape[2+i+j*dx]):
            coords=corner[0]+(i/max_d)*extent, corner[1]-((j+1)/max_d)*extent

            ax.add_patch(Rectangle(coords, extentS, extentS, facecolor=cols[(shape[2+i+j*dx]-1)//4],edgecolor='darkgray',fill=True,hatch=hatchs[(shape[2+i+j*dx]-1)//4],lw=0,transform=ax.transAxes))
            ax.add_patch(Rectangle(coords, extentS,extentS ,edgecolor='gainsboro',fill=False,lw=1.5,transform=ax.transAxes))
            theta=(shape[2+i+j*dx]-1)%4;
            if add_direction:
                ax.arrow(coords[0]+extentS/2, coords[1]+extentS/2, ar_offsets[theta][2]*extentS, ar_offsets[theta][3]*extentS, head_width=0.01*extent, head_length=0.025*extent, fc='k', ec='k',transform=ax.transAxes)

    if rescale:
        ax.text(-0.05, 0.95,f'Size: {sum(i != 0 for i in shape[2:])}\ndX: {dx}\ndY: {dy}',ha='right', va='top',transform = ax.transAxes,fontsize=18)
        #maxB=max(dx,dy)
        #ax.set_xlim([-maxB/2-0.5,maxB/2+0.5])
        #ax.set_ylim([-maxB/2-0.5,maxB/2+0.5])
        ax.set_aspect('equal')
        ax.set_axis_off()
        plt.show(block=False)


faded_lines_alpha=0.3

def plotTransitionsDetailed(pt):
     phen_map_SIZE=defaultdict(dict)
     phen_map_SIZE[1][(1,1,1)]=None
     connection_subsets=defaultdict(list)
     for phen in pt.keys():
          if phen==(1,1,1):
               continue
          for sub_phen in pt[phen].keys():
               connection_subsets[phen].append((phen,sub_phen))

     for k,v in pt.items():
          phen_map_SIZE[np.count_nonzero(k[2:])][k]=v
     counts_map_X={i:len(phen_map_SIZE[v]) for i,v in enumerate(sorted(phen_map_SIZE.keys()))}
     phen_map_SIZE=dict(phen_map_SIZE)
     
     fig,ax = plt.subplots(1)
     connection_dict={}

     phen_dict={(1,1,1):AddPhenotypePatch(ax,(1,1,1),(0,0))}
     
     for i,c in enumerate(sorted(phen_map_SIZE.keys())[1:],1):
          offset=0 if len(phen_map_SIZE[c])%2==1 else .5
          for j,phen in enumerate(sorted(phen_map_SIZE[c])):
               valid_transition=False
               total_weight=sum(phen_map_SIZE[c][phen].values())
                    
               for connector,weight in phen_map_SIZE[c][phen].items():
                    con_size=np.count_nonzero(connector[2:])-1
                    if (con_size+1) not in phen_map_SIZE or connector not in sorted(phen_map_SIZE[con_size+1]):
                         continue
                    valid_transition=True
                    
                    offset2=0 if len(phen_map_SIZE[con_size+1])%2==1 else .5
                    con_y=len(phen_map_SIZE[con_size+1])//2-sorted(phen_map_SIZE[con_size+1]).index(connector)-offset2
                    con_x=sorted(phen_map_SIZE.keys()).index(con_size+1)
                    
                    spline_points=np.array([[con_x+.25,con_y],[con_x+.3,con_y],[i-.3,len(phen_map_SIZE[c])//2-j-offset],[i-.25,len(phen_map_SIZE[c])//2-j-offset]])
                    dx_f=i-con_x
                    rev_tran=dx_f<0
                    dy_f=spline_points[1,1]-spline_points[0,1]
                    dx=(1 if not rev_tran else -1)
                    dy= np.sign(dy_f) if abs(dy_f)>=1 else 0

                    ##special case of "reverse" transition
                    if rev_tran:
                        dx_f=abs(dx_f)
                        spline_points=np.insert(spline_points,2,[con_x,con_y+(0.5 if dy_f>0 else -0.5)],axis=0)
                        dy_f=spline_points[2,1]-spline_points[1,1]
                    steps=1
                    while steps<dx_f:
                         if int(spline_points[steps+(1 if rev_tran else 0),1]*2)%2==counts_map_X[con_x+dx]%2:
                             bump_factor= 0
                         else:
                             bump_factor=.5 if np.sign(spline_points[-1,1]-spline_points[steps+(1 if rev_tran<0 else 0),1])>0 else -.5
                         
                         adjustment_factor=dy+bump_factor
                         if abs(adjustment_factor)>1:
                              adjustment_factor=np.sign(adjustment_factor)*(adjustment_factor%1)
                         spline_points=np.insert(spline_points,steps+(2 if rev_tran else 1),[con_x+dx,spline_points[steps+(1 if rev_tran else 0),1]+adjustment_factor],axis=0)
                         
                         steps+=1
                         dx+=(1 if not rev_tran else -1)
                         dy=dy-(np.sign(dy_f)) if abs(dy)>=1 else 0
                         
                    if rev_tran:     
                         spline_points=np.insert(spline_points,-2,[i,len(phen_map_SIZE[c])//2-j-offset+(.5 if (spline_points[-3,1]-len(phen_map_SIZE[c])//2-j-offset)>0 else-.5)],axis=0)

                    connection_dict[(phen,connector)]=AddConnectionPatch(ax,spline_points,float(weight)/total_weight)

               if valid_transition:
                    phen_dict[phen]=AddPhenotypePatch(ax,phen,(i,len(phen_map_SIZE[c])//2 - j -offset))
        
                    
     ax.set_aspect(1)
     ax.relim()
     ax.autoscale_view()
     ax.grid(False)
     plt.axis('off')

     prev_artists_phens=[]
     prev_artists_lines=[]
     
     def onpick(event):
          coords=event.artist.get_bbox()
          mean_click=np.mean(coords.get_points(),axis=0)
          patch_coord=[int(np.round(mean_click[0])),np.round(mean_click[1]*2)/2.]

          vertical_index=int(counts_map_X[patch_coord[0]]/2-(0 if counts_map_X[patch_coord[0]]%2==1 else .5)-patch_coord[1])
          phen_slices=phen_map_SIZE[sorted(phen_map_SIZE.keys())[patch_coord[0]]]
          phen_key=sorted(phen_slices.keys())[vertical_index]
  
          #! reset colours !#
          for artist in prev_artists_phens:
               artist.set_alpha(.1)
          for artist in prev_artists_lines:
               artist.set_alpha(faded_lines_alpha)
               artist.set_color('gray')
          prev_artists_lines[:] = []
          prev_artists_phens[:] = []

          phen_set={phen_key}
          for phen_pairing, artists in connection_dict.items():
               if phen_key in phen_pairing:
                    phen_set.update(phen_pairing)
                    
                    for artist in artists: #! lines only !#
                         artist.set_alpha(1)
                         artist.set_color('r' if phen_key==phen_pairing[0] else 'b')
                         prev_artists_lines.append(artist)
                         
          for phen in phen_set:
               for artist in phen_dict[phen]:
                    artist.set_alpha(1 if phen==phen_key else 0.4)
                    prev_artists_phens.append(artist)
                    
          fig.canvas.draw()
          return True

     fig.canvas.mpl_connect('pick_event', onpick)
     plt.show(block=False)

##adds a polyomino to a given axis, centered around coordinates xy
def AddPhenotypePatch(ax,shape,xy):
     ar_offsets={0:(0,-.25,0,.25),1:(-.25,0,.25,0),2:(0,.25,0,-.25),3:(.25,0,-.25,0)}
     cols=['darkgreen','royalblue','firebrick','goldenrod','mediumorchid']
     artists=[]
     dx,dy=shape[:2]
     scale=.5/max(shape[:2])
     
     for i,j in product(range(dx),range(dy)):
          if(shape[2+i+j*dx]):
               new_x=xy[0]+(i-dx/2.)*scale
               new_y=xy[1]+(dy/2.-j)*scale-(1.*scale)
               
               artists.append(ax.add_patch(Rectangle((new_x,new_y), scale, scale, fc=cols[(shape[2+i+j*dx]-1)//4],ec='k',fill=True,lw=2,picker=True,alpha=0.1)))
               artists.append(ax.arrow(*(np.array([new_x/scale+.5,new_y/scale+.5,0,0])+ar_offsets[(shape[2+i+j*dx]-1)%4])*scale, head_width=0.075*scale, head_length=0.15*scale, fc='k', ec='k',alpha=0.1))
               
     return artists

def AddConnectionPatch(ax,pts,weight):
     tck, u = splprep(pts.T, u=None, s=0.0,k=3, per=False)
     samples=100
     u_new = np.linspace(u.min(), u.max(), samples)
     x_new, y_new = splev(u_new, tck, der=0)

     ar=ax.arrow(x_new[samples//2],y_new[samples//2],np.diff(x_new[samples//2:samples//2+2])[0],np.diff(y_new[samples//2:samples//2+2])[0], shape='full', lw=0, length_includes_head=True, head_width=.075,alpha=faded_lines_alpha,color='gray')
     ln=ax.plot(x_new, y_new,c='gray', ls='--',lw=weight*2,alpha=faded_lines_alpha)[0]
     
     return (ln,ar)
