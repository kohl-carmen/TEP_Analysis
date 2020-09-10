##------------------------------------------------------------------------------------------------------------
##My HNN core functions
##------------------------------------------------------------------------------------------------------------
## 1) plot_hnn_core_output: Plots output
#    Takes arguments: 
                # dpls:   dipoles (spit our by hnn core)
                # net:    network (spit out by hnn core)
                # name:   figure name (any str)
                # kid:    boolean (display kid's data as well)
                # inputtable:  boolean (display table of input params or not?)
                        # if input_table: further arguments:
                        # params:        dict of all params
                # gaba_table:  boolean (display table of gaba params or not?) (can only select one table)
                        # if gaba_table: further arguments:
                        # param_oi:       list of gaba params
                        # these_params:   values of those gaba params
                        # base_params:    default values to compare these_params to
#    Returns: fig  
                         
## 2) get_dict_param: takes .param file and makes a parameter dictionary
#    Takes arguments: 
                # my_param_out_dir: directory of param file
                # param_name:       name of param file
#    Returns: dict_param  




##------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------
def plot_hnn_core_output(dpls,net,name,kid,input_table,gaba_table,*args):
    if gaba_table==True:
        if input_table==True:
            print('Input Table and GABA Table were both selected. Pick one.')
            gaba_table=False
        try:
            param_oi=args[0]
            these_params=args[1]
            base_params=args[2]
        except:
            print('GABA Table plotting was selected, but no appropriate table arguents were passed. \nGABA Table requires: param_oi, these_params, base_params.')         
            gaba_table=False
    if input_table==True:
        try:
            params=args[0]
        except:
            print('Input Table plotting was selected, but no appropriate table arguents were passed. \nInput Table requires: all_params.')         
            input_table=False      
    
    import numpy as np
    import matplotlib.pyplot as plt
    import pickle
    
    #import data
    data_dir='C:\\Users\\ckohl\\Documents\\SourcePlay\\NewAEP\\'
    data_name= 'adult_right.txt'
    data_dir='C:\\Users\\ckohl\\Documents\\SourcePlay\\NewTiina\\'
    data_name_kid= 'kid(9,10)_right.txt'
    
    
    data_dir='C:\\Users\\ckohl\\Documents\\SourcePlay\\NewTiina\\'
    data_name='Adult_Tone-L Hemi-r.txt'
    data_name_kid='Adult_Tone-R Hemi-r.txt' 
    data_temp=np.loadtxt(data_dir+data_name)
    data_temp_kid=np.loadtxt(data_dir+data_name_kid)
    data_time=np.zeros(len(data_temp))
    data=np.zeros(len(data_temp))
    data_kid=np.zeros(len(data_temp_kid))
    for i in range(0,len(data_temp)):
        data_time[i]=data_temp[i][0]
        data[i]=data_temp[i][1]
        data_kid[i]=data_temp_kid[i][1]
        
    #load original parameters (to compare to)
    dump_dir="C:\\Users\\ckohl\\Desktop\\Current\\HNN Core\\old(wrong data)\\"  
    dump_file='from_here_10_opt_adjusted'
    dump_dir="C:\\Users\\ckohl\\Desktop\\Current\\HNN Core\\"
    dump_file='NEW_manual7_opt10trials'
    dump_file='NEW_manual10_opt_temp5_opt3'
    dump_file='best_new_0114'
    dump_file='Law_best_2401'
    dump_file='Law_best_2401_opt10'
    dump_file='NEW_manual7_opt10trials'
    dump_file='blake_aep_300_scale'
    dump_file='710_tiinascale_opt_pdp'
    f = open(dump_dir+dump_file+'_dpl.txt', 'rb')
    og_dpls= pickle.load(f)
    f.close()
    f = open(dump_dir+dump_file+'_net.txt', 'rb')
    og_net= pickle.load(f)
    f.close()    

    xlim=og_net.params['tstop']
    
    #------------------------------------------------------------------------------------------------------------
    #figure:
    fig = plt.figure(figsize=(20,10))
    #figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()
    #--------------------
    # 1: simulation
    #--------------------
    ax1=fig.add_subplot(221)
    trial=0
    keep_for_mean=np.zeros((len(dpls),len(dpls[0].times)))
    for dpl in dpls:
        ax1.plot(dpl.times,dpl.data['agg'],color=[.5, .5, .5],linewidth=1)
        keep_for_mean[trial,:]=dpl.data['agg']
        trial=trial+1
    ax1.plot(dpl.t,np.mean(keep_for_mean,0),color='k',linewidth=3,label='Model')
    ax1.plot(data_time,data,color='purple',linewidth=3,label='Data Contra')
    if kid==True:
        ax1.plot(data_time,data_kid,color='blue',linewidth=3,label='Data Ipsi')
    ax1.title.set_text('Simulation')
    ax1.set_xlim((0,xlim))
    
    ax1.legend()
        
        
    y_range=ax1.get_ylim()
    y_range=list(y_range)
    
     
    #--------------------
    # 2: spiking    
    #--------------------
    ax2=fig.add_subplot(222)
    #spikes = np.array(sum(net.spiketimes, []))
    #gids = np.array(sum(net.spikegids, []))
    spikes = np.array((net.spikes._times[0]))
    gids = np.array((net.spikes._gids[0]))
    cell_types = ['L5_pyramidal', 'L5_basket', 'L2_pyramidal', 'L2_basket']
    cell_colours=['r','b','g',[.5, .5, .5]]
    count=0
    for cell in cell_types:
        these_spikes=[]
        these_times=[]
        for i in range(0,len(gids)):
            if gids[i]>= net.gid_dict[cell][0] and gids[i]<= net.gid_dict[cell][-1] :
                #if any(gids[i] in s for s in net.gid_dict[cell]):
                these_spikes.append(gids[i])
                these_times.append(spikes[i])
        plt.scatter(these_times,these_spikes,color=cell_colours[count])
        count+=1
    ax2.legend(cell_types)
    ax2.get_yaxis().set_visible(False)
    ax2.title.set_text('Spiking')
    ax2.set_xlim((0,xlim))
    box2=box=ax2.get_position()
    #--------------------
    # 3: dipoles
    #--------------------
    ax3=fig.add_subplot(223)
    trial=0
    keep_for_mean2=np.zeros((len(dpls),len(dpls[0].times)))
    keep_for_mean5=np.zeros((len(dpls),len(dpls[0].times)))
    for dpl in dpls:
        plt.plot(dpl.times,dpl.data['L2'],color='g',linewidth=1)
        plt.plot(dpl.times,dpl.data['L5'],color='r',linewidth=1)
        keep_for_mean2[trial,:]=dpl.data['L2']
        keep_for_mean5[trial,:]=dpl.data['L5']
        trial+=1
    leg1=plt.plot(dpl.times,np.mean(keep_for_mean2,0),color='g',linewidth=3)
    leg1=plt.plot(dpl.times,np.mean(keep_for_mean5,0),color='r',linewidth=3)
    plt.legend(('L2/3','L5'))
    ax3.title.set_text('Dipoles')
    ax3.set_xlim((0,xlim))
    temp_y_range=ax3.get_ylim()
    y_range[0]=min([y_range[0],temp_y_range[0]])
    y_range[1]=max([y_range[1],temp_y_range[1]])
    ax3.set_ylim(y_range)
    #align
    box3=ax3.get_position()
    #--------------------
    # 4: original fit
    #--------------------
    if input_table==False:
        if gaba_table==True:
            ax4=fig.add_subplot(4,4,11)
        else:
            ax4=fig.add_subplot(2,4,7)
        ax4.title.set_text('Original Model')
        trial=0
        keep_for_mean=np.zeros((len(og_dpls),len(og_dpls[0].times)))
        keep_for_mean2=np.zeros((len(og_dpls),len(og_dpls[0].times)))
        keep_for_mean5=np.zeros((len(og_dpls),len(og_dpls[0].times)))
        for dpl in og_dpls:
            keep_for_mean[trial,:]=dpl.data['agg']
            keep_for_mean2[trial,:]=dpl.data['L2']
            keep_for_mean5[trial,:]=dpl.data['L5']
            trial=trial+1
        plt.plot(data_time,data,color='purple',label='data',linewidth=3)
        plt.plot(dpl.times,np.mean(keep_for_mean,0),color='k',label='model',linewidth=3)
        plt.plot(dpl.times,np.mean(keep_for_mean2,0),color='g',linewidth=1)
        plt.plot(dpl.times,np.mean(keep_for_mean5,0),color='r',linewidth=1)
        ax4.set_xlim((0,xlim))
        temp_y_range=ax4.get_ylim()
        y_range[0]=min([y_range[0],temp_y_range[0]])
        y_range[1]=max([y_range[1],temp_y_range[1]])
        #plt.legend(('Data','Model','Layer 2/3','Layer 5'))
        ax4.set_ylim(y_range)
        ax4.get_yaxis().set_visible(False)
        ax3.set_ylim(y_range)
        ax1.set_ylim(y_range)
        #align plots
        box4=ax4.get_position()
        ax4.set_position([box2.x0,box4.y0,box4.width,box4.height])
        if gaba_table==True | input_table==True:
            box4=ax4.get_position()
            ax4.set_position([box4.x0,box3.y0+box3.height-box4.height,box4.width,box4.height])
    #--------------------
    # 5: original spiking
    #--------------------
    if input_table==False:
        if gaba_table==True:
            ax5=fig.add_subplot(4,4,12)
        else:
            ax5=fig.add_subplot(248)
        ax5.title.set_text('Original Spiking')
        #spikes = np.array(sum(og_net.spiketimes, []))
        #gids = np.array(sum(og_net.spikegids, []))
        spikes = np.array((og_net.spikes._times[0]))
        gids = np.array((og_net.spikes._gids[0]))
    
    
        cell_types = ['L5_pyramidal', 'L5_basket', 'L2_pyramidal', 'L2_basket']
        cell_colours=['r','b','g',[.5, .5, .5]]
        count=0
        for cell in cell_types:
            these_spikes=[]
            these_times=[]
            for i in range(0,len(gids)):
                if gids[i]>= og_net.gid_dict[cell][0] and gids[i]<= og_net.gid_dict[cell][-1] :
                    #if any(gids[i] in s for s in net.gid_dict[cell]):
                    these_spikes.append(gids[i])
                    these_times.append(spikes[i])
            plt.scatter(these_times,these_spikes,color=cell_colours[count])
            count+=1
        #ax5.legend(cell_types)
        ax5.get_yaxis().set_visible(False)
        ax5.set_xlim((0,xlim))
        
        #align plots
        box5=ax5.get_position()
        ax5.set_position([box2.x0+box2.width-box5.width,box5.y0,box5.width,box5.height])
        if gaba_table==True:
            box5=ax5.get_position()
            ax5.set_position([box5.x0,box3.y0+box3.height-box5.height,box5.width,box5.height])
        
    if gaba_table==True:
        #--------------------
        # 6: table
        #--------------------
        ax6=fig.add_subplot(4,2,8)
        label_colour=[.8,.8,.8]
        cell_colours=[]
        for p in param_oi:
            cell_colours.append([label_colour,'none','none'])
        collabel=('param','value','base_value')
        rowlabel=[]
        cell_content=[]
        rowlabel=[]
        table_data=np.empty([len(param_oi),2])
        for row in param_oi:
        #make labels
            new_label=[]
            if row[0]=='t':
                new_label='T-'+row[4:8]+row[-1]
            elif row[0]=='s':
                new_label='SD-'+row[10:14]+row[-1]
            else:
                new_label=row[5:]
                new_label=new_label[0:2]+'_I -> '+new_label[9:]
                if new_label[-1]=='b':
                    gaba_type='B'
                elif new_label[-1]=='a':
                    gaba_type='A'
                else:
                    gaba_type='(A)'
                if new_label[10]=='B':
                    new_label=new_label[0:10]+'_I: GABA'+gaba_type
                else:
                    new_label=new_label[0:10]+'_E: GABA'+gaba_type        
            rowlabel.append(new_label)
        #get values
            table_data[len(rowlabel)-1,:]=[these_params[row],base_params[row]]
            cell_content.append([new_label,these_params[row],base_params[row] ])
        #adjust colours
            if these_params[row] < base_params[row]:
                cell_colours[len(rowlabel)-1][1]=[1,.5,.5]
                cell_colours[len(rowlabel)-1][0]=[1,.5,.5]
            elif these_params[row] > base_params[row]:
                cell_colours[len(rowlabel)-1][0]=[.5,1,.5]
                cell_colours[len(rowlabel)-1][1]=[.5,1,.5]

        ax6.axis('off')
        plt.pause(.1)
        table=ax6.table(cellText=cell_content,colLabels=collabel,loc=0,cellColours=cell_colours,cellLoc='center')    
        #table.auto_set_font_size(value=True)
        
        #align plots
        box6=ax6.get_position()
        ax6.set_position([box2.x0,box6.y0,box6.width,box6.height])
        ax6.set_position([box2.x0+.007,box3.y0-.03,box6.width,box6.height])

    if input_table==True:
        #--------------------
        # 6: table
        #--------------------
        ax6=fig.add_subplot(2,2,4)
        ax6.clear()
        # see how many inputs there are
        prox_count=0
        dist_count=0
        collabel=[]
        input_times=[]
        for p in params:
            if p[0:8]=='t_evprox':
                prox_count=prox_count+1
                collabel.append('Prox'+str(prox_count))
                input_times.append(params[p])
            if p[0:8]=='t_evdist':
                dist_count=dist_count+1
                collabel.append('Dist'+str(dist_count))
                input_times.append(params[p])
        #sort inputs
        sort_i=np.argsort(input_times)         
        zipped_pairs = zip(sort_i,collabel)       
        collabel = [x for _, x in sorted(zipped_pairs)]             

        rowlabel=['t','sd','2Pyr_AMPA','2Pyr_NMDA','2Bask_AMPA','2Bas_NMDA','5Pyr_AMPA','5Pyr_NMDA','5Bask_AMPA','5Bask_NMDA']
        actuallabel_pre=['t','sigma_t','gbar','gbar','gbar','gbar','gbar','gbar','gbar','gbar']
        actuallabel_post=['','','_L2Pyr_ampa','_L2Pyr_nmda','_L2Basket_ampa','_L2Basket_nmda','_L5Pyr_ampa','_L5Pyr_nmda','_L5Basket_ampa','_L5Basket_nmda']
        cell_content=[]
        for param_i in range(0,len(rowlabel)):
            this_row=[]
            for input_i in range(0,prox_count+dist_count):
                if collabel[input_i][0]=='P':
                    paramlabel=actuallabel_pre[param_i]+'_evprox_'+collabel[input_i][-1]+actuallabel_post[param_i]
                else:
                    paramlabel=actuallabel_pre[param_i]+'_evdist_'+collabel[input_i][-1]+actuallabel_post[param_i]
                
                try:                  
                    this_row.append(round(params[paramlabel],2))
                except:
                    this_row.append(' ')
                
            this_row.insert(0,rowlabel[param_i])
            cell_content.append(this_row)
            
        collabel.insert(0,' ')
        
        ax6.axis('off')
        plt.pause(.1)
        table=ax6.table(cellText=cell_content,colLabels=collabel,loc=0,cellLoc='center')    
        from matplotlib.font_manager import FontProperties
        for (row, col), cell in table.get_celld().items():
            if (row == 0) or (col == -1):
                cell.set_text_props(fontproperties=FontProperties(weight='bold'))    
        
        #align plots
        box6=ax6.get_position()
        #ax6.set_position([box2.x0+.007,box3.y0-box3.height+.08,box3.width,box6.height*1.8])
        ax6.set_position([box2.x0+.008,box3.y0-.27,box3.width,box3.height*1.8])
    from pylab import suptitle
    fig.suptitle(name, fontsize=14)
    return(fig)








##------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------
 
def get_dict_param(my_param_out_dir,param_name):
    import os.path as op
    dict_param={}
    if op.isfile(my_param_out_dir+'param\\'+param_name +'.param'):#default file in param folder
        convert_this_file=my_param_out_dir+'\\param\\'+param_name +'.param'
    else:
        print(param_name + ' does not exist')    
    #read into dict
    new_f=open(my_param_out_dir +'param_json\\'+param_name +'.json','w+')
    f=open(convert_this_file, 'r')   
    lines = list(f)
    for line in range(0,len(lines)):
        if not (lines[line][0:10]=='sim_prefix' or lines[line][0:10]=='expmt_grou' or lines[line]=='\n' or lines[line]==''):
            div=lines[line].find(':')
            name=lines[line][0:div]
            if name!='spec_cmap':
                value=float(lines[line][div+1:-1])
                new_f.write('"'+name+'":'+ str(value)+',')
                new_f.write('\n')
                dict_param[name]=value
    f.close()
    return(dict_param)
















##------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------
def core_output_basic(dpls,net,name,data,compare,input_table,paramoitable,*args):
    if paramoitable==True:
        if input_table==True:
            print('Input Table and Paramoi Table were both selected. Pick one.')
            paramoitable=False
        try:
            param_oi=args[0] #list of param names
            these_params=net.params
        except:
            print('Paramoi Table plotting was selected, but no appropriate table arguents were passed. \nParamoiable requires: param_oi.')    
            these_params=net.params
            paramoitable=False

    import numpy as np
    import matplotlib.pyplot as plt
    import pickle
    
    #import data
    data_bool=False
    if isinstance(data,str):
        try:
            data_temp=np.loadtxt(data)
            data_bool=True
            data_time=np.zeros(len(data_temp))
            data=np.zeros(len(data_temp))
            for i in range(0,len(data_temp)):
                data_time[i]=data_temp[i][0]
                data[i]=data_temp[i][1]
            
        except:
            print('Data file not found')

    # data_temp2=np.loadtxt('C:\\Users\\ckohl\\Documents\\SourcePlay\\NewTiina\\Adult_Tone-R Hemi-r.txt')
    # data2=np.zeros(len(data_temp2))
    # for i in range(0,len(data_temp2)):
        # data2[i]=data_temp2[i][1]
    
    #load original parameters (to compare to)
    compare_bool=False
    if isinstance(compare,list):
        print('Comparison path given as list (wants string)')
    if isinstance(compare,str):
        try:
            f = open('C:\\Users\\ckohl\\Desktop\\Current\\HNN Core\\'+compare+'_dpl.txt', 'rb')
            og_dpls= pickle.load(f)
            f.close()
            f = open('C:\\Users\\ckohl\\Desktop\\Current\\HNN Core\\'+compare+'_net.txt', 'rb')
            og_net= pickle.load(f)
            f.close()  
            compare_bool=True
            base_params=og_net.params
        except: 
            print('Comparison file not found')
                
     

    if (input_table==True or paramoitable==True) and compare_bool==True:
        print('Parameter tables and Comparisons to previous simulations cannot be plotted simultaneously.')
        input_table=False
        if paramoitable==True:
            compare_bool=False
            compare_for_table=True
    
    if input_table==False and paramoitable==False and compare_bool==False:
        input_table=True
    
    if input_table==True:
        params=net.params
    xlim=net.params['tstop']
    
    #------------------------------------------------------------------------------------------------------------
    #figure:
    fig = plt.figure(figsize=(20,10))
    #figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()
    #--------------------
    # 1: simulation
    #--------------------
    ax1=fig.add_subplot(221)
    trial=0
    keep_for_mean=np.zeros((len(dpls),len(dpls[0].times)))
    for dpl in dpls:
        #ax1.plot(dpl.t,dpl.dpl['agg'],color=[.5, .5, .5],linewidth=1)
        keep_for_mean[trial,:]=dpl.data['agg']
        trial=trial+1
    ax1.plot(dpl.times,np.mean(keep_for_mean,0),color='k',linewidth=3)
    if data_bool==True:
        ax1.plot(data_time,data,color='purple',linewidth=3)
        #ax1.plot(data_time,data2,color='blue',linewidth=3)
    ax1.title.set_text('Simulation')
    ax1.set_xlim((0,xlim))
    ax1.legend(('Model','Data'))
    #ax1.legend(('Model','Contra','Ipsi'))

    y_range=ax1.get_ylim()
    y_range=list(y_range)
    
    #mark inputs
    these_params=net.params
    for temps in these_params:
        if temps[0:2]=='t_':
            if temps[4]=='p':
                Colour='r'
            else:
                Colour='g'
            ax1.arrow(these_params[temps],y_range[0],0,(y_range[1]-y_range[0])/10,color=Colour)
                
            
    #--------------------
    # 2: spiking    
    #--------------------
    ax2=fig.add_subplot(222)
    #spikes = np.array(sum(net.spiketimes, []))
    #gids = np.array(sum(net.spikegids, []))
    spikes = np.array((net.spikes._times[0]))
    gids = np.array((net.spikes._gids[0]))
    cell_types = ['L5_pyramidal', 'L5_basket', 'L2_pyramidal', 'L2_basket']
    cell_colours=['r','b','g',[.5, .5, .5]]
    count=0
    for cell in cell_types:
        these_spikes=[]
        these_times=[]
        for i in range(0,len(gids)):
            if gids[i]>= net.gid_dict[cell][0] and gids[i]<= net.gid_dict[cell][-1] :
                #if any(gids[i] in s for s in net.gid_dict[cell]):
                these_spikes.append(gids[i])
                these_times.append(spikes[i])
        plt.scatter(these_times,these_spikes,color=cell_colours[count])
        count+=1
    ax2.legend(cell_types)
    ax2.get_yaxis().set_visible(False)
    ax2.title.set_text('Spiking')
    ax2.set_xlim((0,xlim))
    box2=box=ax2.get_position()
    #--------------------
    # 3: dipoles
    #--------------------
    ax3=fig.add_subplot(223)
    trial=0
    keep_for_mean2=np.zeros((len(dpls),len(dpls[0].times)))
    keep_for_mean5=np.zeros((len(dpls),len(dpls[0].times)))
    for dpl in dpls:
        plt.plot(dpl.times,dpl.data['L2'],color='g',linewidth=1)
        plt.plot(dpl.times,dpl.data['L5'],color='r',linewidth=1)
        keep_for_mean2[trial,:]=dpl.data['L2']
        keep_for_mean5[trial,:]=dpl.data['L5']
        trial+=1
    leg1=plt.plot(dpl.times,np.mean(keep_for_mean2,0),color='g',linewidth=3)
    leg1=plt.plot(dpl.times,np.mean(keep_for_mean5,0),color='r',linewidth=3)
    plt.legend(('L2/3','L5'))
    ax3.title.set_text('Dipoles')
    ax3.set_xlim((0,xlim))
    temp_y_range=ax3.get_ylim()
    y_range[0]=min([y_range[0],temp_y_range[0]])
    y_range[1]=max([y_range[1],temp_y_range[1]])
    ax3.set_ylim(y_range)
    #align
    box3=ax3.get_position()
    #--------------------
    # 4: original fit
    #--------------------
    if compare_bool==True:
        ax4=fig.add_subplot(2,4,7)
        ax4.title.set_text('Original Model')
        trial=0
        keep_for_mean=np.zeros((len(og_dpls),len(og_dpls[0].times)))
        keep_for_mean2=np.zeros((len(og_dpls),len(og_dpls[0].times)))
        keep_for_mean5=np.zeros((len(og_dpls),len(og_dpls[0].times)))
        for dpl in og_dpls:
            keep_for_mean[trial,:]=dpl.data['agg']
            keep_for_mean2[trial,:]=dpl.data['L2']
            keep_for_mean5[trial,:]=dpl.data['L5']
            trial=trial+1
            if data_bool==True:
                plt.plot(data_time,data,color='purple',label='data',linewidth=3)
        plt.plot(dpl.times,np.mean(keep_for_mean,0),color='k',label='model',linewidth=3)
        plt.plot(dpl.times,np.mean(keep_for_mean2,0),color='g',linewidth=1)
        plt.plot(dpl.times,np.mean(keep_for_mean5,0),color='r',linewidth=1)
        ax4.set_xlim((0,xlim))
        temp_y_range=ax4.get_ylim()
        y_range[0]=min([y_range[0],temp_y_range[0]])
        y_range[1]=max([y_range[1],temp_y_range[1]])
        #plt.legend(('Data','Model','Layer 2/3','Layer 5'))
        ax4.set_ylim(y_range)
        ax4.get_yaxis().set_visible(False)
        ax3.set_ylim(y_range)
        ax1.set_ylim(y_range)
        #align plots
        box4=ax4.get_position()
        ax4.set_position([box2.x0,box4.y0,box4.width,box4.height])
        if paramoitable==True | input_table==True:
            box4=ax4.get_position()
            ax4.set_position([box4.x0,box3.y0+box3.height-box4.height,box4.width,box4.height])
    #--------------------
    # 5: original spiking
    #--------------------
    if compare_bool==True:
        ax5=fig.add_subplot(2,4,8)
        ax5.title.set_text('Original Spiking')
        spikes = np.array(sum(og_net.spikes._times, []))
        gids = np.array(sum(og_net.spikes._gids, []))
        cell_types = ['L5_pyramidal', 'L5_basket', 'L2_pyramidal', 'L2_basket']
        cell_colours=['r','b','g',[.5, .5, .5]]
        count=0
        for cell in cell_types:
            these_spikes=[]
            these_times=[]
            for i in range(0,len(gids)):
                if gids[i]>= og_net.gid_dict[cell][0] and gids[i]<= og_net.gid_dict[cell][-1] :
                    #if any(gids[i] in s for s in net.gid_dict[cell]):
                    these_spikes.append(gids[i])
                    these_times.append(spikes[i])
            plt.scatter(these_times,these_spikes,color=cell_colours[count])
            count+=1
        #ax5.legend(cell_types)
        ax5.get_yaxis().set_visible(False)
        ax5.set_xlim((0,xlim))
        
        #align plots
        box5=ax5.get_position()
        ax5.set_position([box2.x0+box2.width-box5.width,box5.y0,box5.width,box5.height])
        if paramoitable==True:
            box5=ax5.get_position()
            ax5.set_position([box5.x0,box3.y0+box3.height-box5.height,box5.width,box5.height])
        
    if paramoitable==True:
        #--------------------
        # 6: table
        #--------------------
        ax4=fig.add_subplot(2,2,4)
        label_colour=[.8,.8,.8]
        cell_colours=[]
        for p in param_oi:
            cell_colours.append([label_colour,'none','none'])
        if compare_for_table==True:
            collabel=('param','value','base_value')
        else:
            collabel=('param','value')
        rowlabel=[]
        cell_content=[]
        rowlabel=[]
        table_data=np.empty([len(param_oi),2])
        for row in param_oi:    
            rowlabel.append(row)
            if compare_for_table==True:
                cell_content.append([row,these_params[row],base_params[row] ])
                
                if these_params[row] < base_params[row]:
                    cell_colours[len(rowlabel)-1][1]=[1,.5,.5]
                    cell_colours[len(rowlabel)-1][0]=[1,.5,.5]
                elif these_params[row] > base_params[row]:
                    cell_colours[len(rowlabel)-1][0]=[.5,1,.5]
                    cell_colours[len(rowlabel)-1][1]=[.5,1,.5]
            else:
                cell_content.append([row,these_params[row]])

                
        #adjust colours
        if compare_for_table==True:
            ax4.axis('off')
            plt.pause(.1)
            table=ax4.table(cellText=cell_content,colLabels=collabel,loc=0,cellColours=cell_colours,cellLoc='center')    
            #table.auto_set_font_size(value=True)
        else:
            ax4.axis('off')
            plt.pause(.1)
            table=ax4.table(cellText=cell_content,colLabels=collabel,loc=0,cellLoc='center')    
       
        
        #align plots
        box4=ax4.get_position()
        ax4.set_position([box2.x0+.007,box3.y0-.03,box4.width,box4.height])

    if input_table==True:
        #--------------------
        # 6: table
        #--------------------
        ax4=fig.add_subplot(2,2,4)
        ax4.clear()
        # see how many inputs there are
        prox_count=0
        dist_count=0
        collabel=[]
        input_times=[]
        for p in params:
            if p[0:8]=='t_evprox':
                prox_count=prox_count+1
                collabel.append('Prox'+str(prox_count))
                input_times.append(params[p])
            if p[0:8]=='t_evdist':
                dist_count=dist_count+1
                collabel.append('Dist'+str(dist_count))
                input_times.append(params[p])
        #sort inputs
        sort_i=np.argsort(input_times)         
        zipped_pairs = zip(sort_i,collabel)       
        collabel = [x for _, x in sorted(zipped_pairs)]             

        rowlabel=['t','sd','2Pyr_AMPA','2Pyr_NMDA','2Bask_AMPA','2Bas_NMDA','5Pyr_AMPA','5Pyr_NMDA','5Bask_AMPA','5Bask_NMDA']
        actuallabel_pre=['t','sigma_t','gbar','gbar','gbar','gbar','gbar','gbar','gbar','gbar']
        actuallabel_post=['','','_L2Pyr_ampa','_L2Pyr_nmda','_L2Basket_ampa','_L2Basket_nmda','_L5Pyr_ampa','_L5Pyr_nmda','_L5Basket_ampa','_L5Basket_nmda']
        cell_content=[]
        for param_i in range(0,len(rowlabel)):
            this_row=[]
            for input_i in range(0,prox_count+dist_count):
                if collabel[input_i][0]=='P':
                    paramlabel=actuallabel_pre[param_i]+'_evprox_'+collabel[input_i][-1]+actuallabel_post[param_i]
                else:
                    paramlabel=actuallabel_pre[param_i]+'_evdist_'+collabel[input_i][-1]+actuallabel_post[param_i]
                
                try:                  
                    this_row.append(round(params[paramlabel],2))
                except:
                    this_row.append(' ')
                
            this_row.insert(0,rowlabel[param_i])
            cell_content.append(this_row)
            
        collabel.insert(0,' ')
        
        ax4.axis('off')
        plt.pause(.1)
        table=ax4.table(cellText=cell_content,colLabels=collabel,loc=0,cellLoc='center')    
        from matplotlib.font_manager import FontProperties
        for (row, col), cell in table.get_celld().items():
            if (row == 0) or (col == -1):
                cell.set_text_props(fontproperties=FontProperties(weight='bold'))    
        
        #align plots
        box6=ax4.get_position()
        #ax6.set_position([box2.x0+.007,box3.y0-box3.height+.08,box3.width,box6.height*1.8])
        ax4.set_position([box2.x0+.008,box3.y0-.27,box3.width,box3.height*1.8])
    from pylab import suptitle
    fig.suptitle(name, fontsize=14)
    plt.pause(.1)
    return(fig)