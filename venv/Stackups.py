import pandas as pd
import numpy as np

pd.set_option('display.expand_frame_repr', False)

k = 1.3806485279e-23
T0 = 290

class RF_stackup():

    def __init__(self,names = 'thing1',T0 = 290):
        self.stack = pd.DataFrame(columns=['Type','Gain_dB', 'NF_dB', 'ITOI_dB', 'IP2_dB', 'OTOI_dB','IP1dB', 'DANL', 'VSWR'])
        self.stack.ix[:,'Type'] = names
        self.T0 = T0
        self.flags = 0x0000 #this number is the initialized flag
        self.Gain_dB_mask = 0x0001
        self.NF_dB_mask = 0x0010
        self.Gain_Lin_mask = 0x0002
        self.Cum_Gain_dB_mask = 0x0004
        self.Cum_Gain_Lin_mask = 0x0008
        self.NF_Lin_mask = 0x0020
        self.cascade_NF_mask = 0x0040
        self.ITOI_dB_mask = 0x0080
        self.OTOI_lin_mask = 0x0100
        self.OTOI_dB_mask = 0x0200
        self.Coh_Gain_mask = 0x0400
        self.VSWR_mask = 0x0800
        self.IP2_dB_mask = 0x1000
        self.OP2_lin_mask = 0x2000

    def populate_gains_dB(self,gains_dB):

        space = len(self.stack['Gain_dB'])
        values = len(gains_dB)

        if space == values:
            self.stack['Gain_dB'] = gains_dB
        else:
            raise ValueError('Error(populate gains): there are ', space, ' spaces in the array and you gave me ', values, ' values')

        if bin(self.flags & self.Gain_dB_mask).count("1") == 0:
            self.flags = self.flags + self.Gain_dB_mask


        self.calculate_linear_gain()


    def calculate_linear_gain(self):
        self.stack["Gain_Lin"] = 10**(self.stack['Gain_dB']/10)

        if bin(self.flags & self.Gain_Lin_mask).count("1") == 0:
            self.flags = self.flags + self.Gain_Lin_mask

    def calculate_Cum_gain_dB(self):

        self.calculate_art()

        temp = self.stack.copy()
        temp_gain = temp['Gain_dB']

        for i in range(0,len(self.stack["Type"])):

            if self.stack.ix[i,"Type"] == 'combiner':
                temp.ix[i,'Gain_dB'] = self.stack.ix[i,'Gain_dB'] # + 10*np.log10(self.stack.ix[i,'Device Inputs'])


        self.stack["Cum_Gain_dB"] = np.cumsum(temp['Gain_dB'])

        #self.stack['Min_Gain_dB'] = np.zeros_like(self.stack['Gain_dB'])
        temp = self.stack.copy()
        temp_gain = temp['Gain_dB']
        #
        for i in range(0,len(self.stack["Type"])):
             if self.stack.ix[i,"Type"] == 'combiner':
                 temp.ix[i,'Gain_dB'] = self.stack.ix[i,'Gain_dB'] #+ 10*np.log10(self.stack.ix[i,'Device Inputs'])
             if self.stack.ix[i,"Type"] == 'cable' or self.stack.ix[i,"Type"] == 'tx_line':
                 temp.ix[i,'Gain_dB'] = self.stack.ix[i,'|acbl|min_dB'] #+ 10*np.log10(self.stack.ix[i,'Device Inputs'])
        #
        self.stack["Cum_Min_Gain_dB"] = np.cumsum(temp['Gain_dB'])

        #self.stack['Max_Gain_dB'] = np.zeros_like(self.stack['Gain_dB'])
        temp = self.stack.copy()
        temp_gain = temp['Gain_dB']
        #
        for i in range(0, len(self.stack["Type"])):
            if self.stack.ix[i, "Type"] == 'combiner':
                temp.ix[i, 'Gain_dB'] = self.stack.ix[i, 'Gain_dB'] #+ 10 * np.log10(self.stack.ix[i, 'Device Inputs'])
            if self.stack.ix[i, "Type"] == 'cable' or self.stack.ix[i, "Type"] == 'tx_line':
                temp.ix[i, 'Gain_dB'] = self.stack.ix[i, '|acbl|max_dB'] #+ 10 * np.log10(self.stack.ix[i, 'Device Inputs'])
        #
        self.stack["Cum_Max_Gain_dB"] = np.cumsum(temp['Gain_dB'])



        if bin(self.flags & self.Cum_Gain_dB_mask).count("1") == 0:
            self.flags = self.flags + self.Cum_Gain_dB_mask

    def calculate_Cum_gain_linear(self):

        if (self.flags & self.Gain_dB_mask) == 1:

            temp = self.stack.copy()
            temp_gain = temp['Gain_dB']

            for i in range(0, len(self.stack["Type"])):

                if self.stack.ix[i, "Type"] == 'combiner':
                    temp.ix[i,'Gain_dB'] = self.stack.ix[i, 'Gain_dB'] + 10*np.log10(2)

            self.stack["Cum_Gain_Lin"] = 10**(np.cumsum(temp_gain)/10)

        elif (self.flags & self.Gain_Lin_mask) == 1:

            temp = self.stack.copy()
            temp_gain = temp['Cum_Gain_Lin']

            for i in range(0, len(self.stack["Type"])):

                if self.stack.ix[i, "Type"] == 'combiner':
                    temp_gain[i] = 2*temp_gain

            self.stack["Cum_Gain_lin"] = np.cumprod(temp_gain)


        else:
            raise ValueError('Error: No Gain Values to work with')

        if bin(self.flags & self.Cum_Gain_Lin_mask).count("1") == 0:
            self.flags = self.flags + self.Cum_Gain_Lin_mask


    def populate_NF_dB(self,NF_dB):

        space = len(self.stack['NF_dB'])
        values = len(NF_dB)

        if space == values:
            self.stack['NF_dB'] = NF_dB
            if any(n<0 for n in NF_dB):
                raise ValueError('Error: Negative number in NF_dB, this is not possible!!')
        else:
            raise ValueError('Error(populate NF): there are ', space, ' spaces in the array and you gave me ', values, ' values')

        if bin(self.flags & self.NF_dB_mask).count("1") == 0:
            self.flags = self.flags + self.NF_dB_mask

    def calculate_linear_NF(self):


        if bin(self.flags & self.NF_dB_mask).count("1") == 1:
            self.stack["NF_Lin"] = 10**(self.stack['NF_dB']/10)


        else:
            raise ValueError('Error: Need to populate NF values in dB')

        if bin(self.flags & self.NF_Lin_mask).count("1") == 0:
            self.flags = self.flags + self.NF_Lin_mask

    def cascade_NF(self):

        if bin(self.flags & self.Cum_Gain_Lin_mask).count("1") == 0:
            if bin(self.flags & self.Cum_Gain_dB_mask).count("1") == 0:
                self.calculate_Cum_gain_dB()
                self.calculate_Cum_gain_linear()
            else:
                self.calculate_Cum_gain_linear()

        if bin(self.flags & self.NF_Lin_mask).count("1") == 0:
            self.calculate_linear_NF()


        zero_ones = np.ones_like(self.stack["NF_Lin"])
        zero_ones[0] = 0

        adj_gains = np.roll(self.stack["Cum_Gain_Lin"],1)
        adj_gains[0] = 1

        fbs = (self.stack["NF_Lin"] - zero_ones)/adj_gains

        nf = np.cumsum(fbs)

        self.stack["Cum_NF_Lin"] = nf
        self.stack["Cum_NF_dB"] = 10*np.log10(nf)

        mixer_locations = self.stack.Type[self.stack.Type == 'mixer'].index.tolist()
        filter_locations = self.stack.Type[self.stack.Type == 'filter'].index.tolist()
        mixer_ssb = self.stack.Type[self.stack.Type == 'mixer_wb'].index.tolist()
        mixer_dsb = self.stack.Type[self.stack.Type == 'mixer_sup'].index.tolist()


        fb_gp_prime = np.zeros_like(self.stack["NF_Lin"])



        if len(mixer_locations) > 0:
            for i in range(0, len(mixer_locations)):
                if len(filter_locations) > 0:
                    filters_in_range = [x for x in filter_locations if x<mixer_locations[i]]

                    if len(filters_in_range)==0:
                        fb_gp_prime[mixer_locations[i]] = (1 + (
                                    self.stack.ix[mixer_locations[i] - 1, 'Cum_NF_Lin'] - 1)) * (self.stack.ix[
                                                              mixer_locations[i] - 1, 'Cum_Gain_Lin']-1)

                    else:
                        filt_ind = max(filters_in_range)

                        fb_gp_prime[mixer_locations[i]] = (1 + (self.stack.ix[mixer_locations[i] - 1, 'Cum_NF_Lin'] - self.stack.ix[filt_ind, 'Cum_NF_Lin'])*self.stack.ix[filt_ind, 'Cum_Gain_Lin']) * ((self.stack.ix[
                        mixer_locations[i] - 1, 'Cum_Gain_Lin']-1)/self.stack.ix[filt_ind, 'Cum_Gain_Lin'])





                else:
                    fb_gp_prime[mixer_locations[i]] = (1 + (self.stack.ix[mixer_locations[i] - 1, 'Cum_NF_Lin'] - 1)) * (self.stack.ix[
                        mixer_locations[i] - 1, 'Cum_Gain_Lin']-1)

        nf_eq = self.stack["NF_Lin"] + fb_gp_prime



        if len(mixer_ssb) > 0:
            for i in range(0,len(mixer_ssb)):
                nf_eq[mixer_ssb[i]] = 2 * (self.stack.ix[mixer_ssb[i], 'NF_Lin']-1) + 1 + 1

        if len(mixer_dsb) > 0:
            for i in range(0,len(mixer_dsb)):
                nf_eq[mixer_dsb[i]] = 2 * (self.stack.ix[mixer_dsb[i], 'NF_Lin']-1) + 1 + 0

        self.stack["Eq_NF_Lin"] = nf_eq
        self.stack["Eq_NF_dB"] = 10 * np.log10(nf_eq)

        zero_ones = np.ones_like(self.stack["Eq_NF_Lin"])
        zero_ones[0] = 0

        adj_gains = np.roll(self.stack["Cum_Gain_Lin"], 1)
        adj_gains[0] = 1

        fbs_eq = (self.stack["Eq_NF_Lin"] - zero_ones) / adj_gains

        nf_eq = np.cumsum(fbs_eq)

        self.stack["Eq_Cum_NF_Lin"] = nf_eq
        self.stack["Eq_Cum_NF_dB"] = 10 * np.log10(nf_eq)

        if bin(self.flags & self.cascade_NF_mask).count("1") == 0:
            self.flags = self.flags + self.cascade_NF_mask


    def populate_ITOI_dB(self,ITOI_dB):

        space = len(self.stack['ITOI_dB'])
        values = len(ITOI_dB)

        if space == values:
            self.stack['ITOI_dB'] = ITOI_dB
        else:
            raise ValueError('Error(populate gains): there are ', space, ' spaces in the array and you gave me ', values, ' values')

        if bin(self.flags & self.ITOI_dB_mask).count("1") == 0:
            self.flags = self.flags + self.ITOI_dB_mask

    def calculate_OTOI_dB(self):
        if bin(self.flags & self.Gain_dB_mask).count("1") == 1:
            if bin(self.flags & self.ITOI_dB_mask).count("1") == 1:
                self.stack.ix[:,'ITOI_dB'].fillna(1000,inplace=True)
                self.stack['OTOI_dB'] = self.stack['ITOI_dB'] + self.stack['Gain_dB']
        else:
            raise ValueError('Error: cannot populate OTOI if there is no Gain available')
        self.calculate_lin_OTOI()

    def calculate_lin_OTOI(self):
        self.stack['OTOI_Lin'] = 10**(self.stack['OTOI_dB']/10)

        if bin(self.flags & self.OTOI_lin_mask).count("1") == 0:
            self.flags = self.flags + self.OTOI_lin_mask

    def cascade_OTOI(self):
        if bin(self.flags & self.OTOI_lin_mask).count("1") == 1:

            if bin(self.flags & self.Gain_Lin_mask).count("1") == 0:
                self.calculate_linear_gain()

            self.stack["Cum_OTOI_Lin"] = np.zeros_like(self.stack["OTOI_Lin"])

            gain = self.stack["Gain_Lin"].values
            elem_OTOI = self.stack["OTOI_Lin"].values
            old_otoi = 0

            for x in range(0,len(self.stack["OTOI_Lin"])):
                if x == 0:
                    self.stack.ix[0,"Cum_OTOI_Lin"] = elem_OTOI[x]
                    old_otoi = elem_OTOI[x]
                else:
                    self.stack.ix[x,"Cum_OTOI_Lin"] = 1/((1/elem_OTOI[x])+(1/(old_otoi*gain[x])))
                    old_otoi = 1/((1/elem_OTOI[x])+(1/(old_otoi*gain[x])))

            self.stack["Cum_OTOI_dB"] = 10 * np.log10(self.stack["Cum_OTOI_Lin"])
            self.stack["Cum_ITOI_dB"] = 10 * np.log10(self.stack["Cum_OTOI_Lin"]) - self.stack["Cum_Gain_dB"]

        else:
            raise ValueError('Error: cannot cascade linear OTOI if there is no linear OTOI data available')

        if bin(self.flags & self.OTOI_dB_mask).count("1") == 0:
            self.flags = self.flags + self.OTOI_dB_mask

    def populate_IP2_dB(self,IP2_dB):

        space = len(self.stack['IP2_dB'])
        values = len(IP2_dB)

        if space == values:
            self.stack['IP2_dB'] = IP2_dB
        else:
            raise ValueError('Error(populate ip2): there are ', space, ' spaces in the array and you gave me ', values, ' values')

        if bin(self.flags & self.IP2_dB_mask).count("1") == 0:
            self.flags = self.flags + self.IP2_dB_mask

    def calculate_OP2_dB(self):
        if bin(self.flags & self.Gain_dB_mask).count("1") == 1:
            if bin(self.flags & self.IP2_dB_mask).count("1") == 1:
                self.stack.ix[:,'IP2_dB'].fillna(1000,inplace=True)
                self.stack['OP2_dB'] = self.stack['IP2_dB'] + self.stack['Gain_dB']
        else:
            raise ValueError('Error: cannot populate OP2 if there is no Gain available')
        self.calculate_lin_OP2()

    def calculate_lin_OP2(self):
        self.stack['OP2_Lin'] = 10**(self.stack['OP2_dB']/10)

        if bin(self.flags & self.OP2_lin_mask).count("1") == 0:
            self.flags = self.flags + self.OP2_lin_mask

    def cascade_OP2(self):
        if bin(self.flags & self.OP2_lin_mask).count("1") == 1:

            if bin(self.flags & self.Gain_Lin_mask).count("1") == 0:
                self.calculate_linear_gain()

            self.stack["Cum_OP2_Lin"] = np.zeros_like(self.stack["OP2_Lin"])

            gain = self.stack["Gain_Lin"].values
            elem_OP2 = self.stack["OP2_Lin"].values
            old_op2 = 0

            for x in range(0,len(self.stack["OP2_Lin"])):
                if x == 0:
                    self.stack.ix[0,"Cum_OP2_Lin"] = elem_OP2[x]
                    old_op2 = elem_OP2[x]
                else:
                    self.stack.ix[x,"Cum_OP2_Lin"] = 1/((np.sqrt(1/elem_OP2[x])+np.sqrt(1/(old_op2*gain[x])))**2)
                    old_op2 = 1/((np.sqrt(1/elem_OP2[x])+np.sqrt(1/(old_op2*gain[x])))**2)

            self.stack["Cum_OP2_dB"] = 10 * np.log10(self.stack["Cum_OP2_Lin"])
            self.stack["Cum_IP2_dB"] = 10 * np.log10(self.stack["Cum_OP2_Lin"]) - self.stack["Cum_Gain_dB"]

        else:
            raise ValueError('Error: cannot cascade linear OTOI if there is no linear OTOI data available')

        if bin(self.flags & self.OTOI_dB_mask).count("1") == 0:
            self.flags = self.flags + self.OTOI_dB_mask

    def calculate_OP1dB_dB(self):
        self.stack.ix[:,'IP1dB'].fillna(1000,inplace=True)
        print(self.stack['Gain_dB'])
        self.stack['OP1dB'] = self.stack['IP1dB'] + self.stack['Gain_dB']
        self.calculate_OP1dB_lin()

    def calculate_OP1dB_lin(self):
        self.stack['OP1dB_Lin'] = 10 ** (self.stack['OP1dB'] / 10)
        self.cascade_OP1dB_lin()

    def cascade_OP1dB_lin(self):
        if bin(self.flags & self.Gain_Lin_mask).count("1") == 0:
            self.calculate_linear_gain()
        if bin(self.flags & self.Cum_Gain_Lin_mask).count("1") == 0:
            self.calculate_Cum_gain_linear()
        if bin(self.flags & self.Cum_Gain_dB_mask).count("1") == 0:
            self.calculate_Cum_gain_dB()

        self.stack["Cum_OP1dB_Lin"] = np.zeros_like(self.stack["OP1dB_Lin"])

        gain = self.stack["Gain_Lin"].values
        elem_OP1dB = self.stack["OP1dB_Lin"].values
        old_op1db = 0

        for x in range(0,len(self.stack["OP1dB_Lin"])):
            if x == 0:
                self.stack.ix[0,"Cum_OP1dB_Lin"] = elem_OP1dB[0]
                old_op1db = elem_OP1dB[0]
            else:
                self.stack.ix[x,"Cum_OP1dB_Lin"] = 1/((1/elem_OP1dB[x])+(1/(old_op1db*gain[x])))
                old_op1db = 1/((1/elem_OP1dB[x])+(1/(old_op1db*gain[x])))

        self.stack["Cum_OP1dB_dB"] = 10 * np.log10(self.stack["Cum_OP1dB_Lin"])
        self.stack["Cum_IP1dB_dB"] = 10 * np.log10(self.stack["Cum_OP1dB_Lin"]) - (self.stack["Cum_Gain_dB"])

    def calculate_Coherent_Gain(self,inputs):

        self.stack["Coh_Gain_dB"] = self.stack['Gain_dB'] + 20*np.log10(inputs)

        ####################################
        # calculate min coherent gain
        #####################################

        temp = self.stack.copy()
        #
        for i in range(0, len(self.stack["Type"])):
            if self.stack.ix[i, "Type"] == 'combiner':
                temp.ix[i, 'Gain_dB'] = self.stack.ix[i, 'Gain_dB'] + 10 * np.log10(self.stack.ix[i, 'Device Inputs'])
            if self.stack.ix[i, "Type"] == 'cable' or self.stack.ix[i, "Type"] == 'tx_line':
                temp.ix[i, 'Gain_dB'] = self.stack.ix[i, '|acbl|min_dB'] + 10 * np.log10(
                    self.stack.ix[i, 'Device Inputs'])
        #

        self.stack["Min_Coh_Gain_dB"] = temp['Gain_dB'] + 20 * np.log10(inputs)

        ####################################
        # calculate max coherent gain
        #####################################

        temp = self.stack.copy()
        #
        for i in range(0, len(self.stack["Type"])):
            if self.stack.ix[i, "Type"] == 'combiner':
                temp.ix[i, 'Gain_dB'] = self.stack.ix[i, 'Gain_dB']
            if self.stack.ix[i, "Type"] == 'cable' or self.stack.ix[i, "Type"] == 'tx_line':
                temp.ix[i, 'Gain_dB'] = self.stack.ix[i, '|acbl|max_dB']
        #

        self.stack["Max_Coh_Gain_dB"] = temp['Gain_dB'] + 20 * np.log10(inputs)

        self.stack["Device Inputs"] = inputs

        for i in range(0,len(self.stack["Coh_Gain_dB"])):

            if self.stack.ix[i,"Type"] == 'combiner':
                self.stack.ix[i,"Coh_Gain_dB"] = self.stack.ix[i,"Coh_Gain_dB"]
                self.stack.ix[i, "Min_Coh_Gain_dB"] = self.stack.ix[i, "Min_Coh_Gain_dB"]
                self.stack.ix[i, "Max_Coh_Gain_dB"] = self.stack.ix[i, "Max_Coh_Gain_dB"]
            else:
                self.stack.ix[i, "Coh_Gain_dB"] = self.stack.ix[i, "Coh_Gain_dB"]
                self.stack.ix[i, "Min_Coh_Gain_dB"] = self.stack.ix[i, "Min_Coh_Gain_dB"]
                self.stack.ix[i, "Max_Coh_Gain_dB"] = self.stack.ix[i, "Max_Coh_Gain_dB"] 


        self.stack["Cum_Coh_Gain_dB"] = np.cumsum(self.stack["Coh_Gain_dB"]) #+ 20 * np.log10(inputs)

        if bin(self.flags & self.Coh_Gain_mask).count("1") == 0:
            self.flags = self.flags + self.Coh_Gain_mask

    def calculate_siginout(self,input_power_dBm):

        if bin(self.flags & self.Gain_dB_mask).count("1") == 1:

            #self.stack["Signal_Out_dBm"] = np.zeros_like(self.stack["Gain_dB"])

            temp = self.stack.copy()

            temp = temp["Coh_Gain_dB"]

            temp[0] = temp[0] + input_power_dBm

            temp = np.cumsum(temp)

            self.stack["Signal_In_dBm"] = np.zeros_like(temp)
            self.stack["Signal_Out_dBm"] = np.zeros_like(temp)
            self.stack["Compressing"] = np.zeros_like(temp)
            self.stack["Actual_Gain_dB"] = np.zeros_like(temp)

            old_gain = input_power_dBm

            for i in range(0,len(temp)):

                self.stack.ix[i, "Signal_In_dBm"] = old_gain

                if self.stack.ix[i, "Signal_In_dBm"]+self.stack.ix[i, "Gain_dB"] > self.stack.ix[i,"IP1dB"]:
                    self.stack.ix[i, "Signal_Out_dBm"]  = self.stack.ix[i,"IP1dB"]
                    self.stack.ix[i,"Compressing"] = 1
                elif self.stack.ix[i, "Signal_In_dBm"]+self.stack.ix[i, "Gain_dB"] > self.stack.ix[i,"IP1dB"]-3:
                    self.stack.ix[i, "Signal_Out_dBm"] = self.stack.ix[i, "Signal_In_dBm"]+self.stack.ix[i, "Gain_dB"]
                    self.stack.ix[i,"Compressing"] = 0.5
                else:
                    self.stack.ix[i, "Signal_Out_dBm"] = self.stack.ix[i, "Signal_In_dBm"] + self.stack.ix[i, "Gain_dB"]

                old_gain = self.stack.ix[i, "Signal_Out_dBm"]


            self.stack["Actual_Gain_dB"] =  self.stack["Signal_Out_dBm"]-self.stack["Signal_In_dBm"]
            self.stack["Cum_Actual_Gain_dB"] = np.cumsum(self.stack["Actual_Gain_dB"])
            self.stack["IP1dB_Headroom_dB"] = self.stack["IP1dB"]-self.stack["Signal_Out_dBm"]

        else:

            raise ValueError('Error: Gain parameters must be loaded before sigout calc')

    def calculate_noise_gain(self):

        self.stack["Noise_Gain_dB"]  = self.stack["Actual_Gain_dB"] + 10*np.log10(self.stack["Device Inputs"])
        self.stack["Noise_Gain_Lin"] = 10**(self.stack["Noise_Gain_dB"]/10)
        self.calculate_lin_noise_added()

    def calculate_lin_noise_added(self):
        self.stack["Eq_Noise_added_Lin"] = k*T0*(self.stack["Eq_NF_Lin"]-1)*self.stack["Noise_Gain_Lin"]
        self.stack["Noise_added_Lin"] = k * T0 * (self.stack["NF_Lin"] - 1) * self.stack["Noise_Gain_Lin"]

        temp = self.stack.copy()
        temp = temp["Eq_Noise_added_Lin"]
        self.stack["Eq_Noise_Out_Lin_Hz"] = np.zeros_like(temp)
        self.stack["Eq_Noise_In_Lin_Hz"] = np.zeros_like(temp)
        self.stack["Noise_Out_Lin_Hz"] = np.zeros_like(temp)
        self.stack["Noise_In_Lin_Hz"] = np.zeros_like(temp)

        for i in range(0, len(temp)):
            if i == 0:
                self.stack.ix[i,"Eq_Noise_Out_Lin_Hz"] = k * T0 * self.stack.ix[i,"Noise_Gain_Lin"] + self.stack.ix[i,"Eq_Noise_added_Lin"]
                self.stack.ix[i, "Eq_Noise_In_Lin_Hz"] = k * T0

                self.stack.ix[i, "Noise_Out_Lin_Hz"] = k * T0 * self.stack.ix[i, "Noise_Gain_Lin"] + self.stack.ix[
                    i, "Noise_added_Lin"]
                self.stack.ix[i, "Noise_In_Lin_Hz"] = k * T0

            else:
                self.stack.ix[i, "Eq_Noise_Out_Lin_Hz"] = self.stack.ix[i-1,"Eq_Noise_Out_Lin_Hz"] * self.stack.ix[i, "Noise_Gain_Lin"] + self.stack.ix[i,"Eq_Noise_added_Lin"]
                self.stack.ix[i, "Eq_Noise_In_Lin_Hz"]  = self.stack.ix[i-1,"Eq_Noise_Out_Lin_Hz"]

                self.stack.ix[i, "Noise_Out_Lin_Hz"] = self.stack.ix[i - 1, "Noise_Out_Lin_Hz"] * self.stack.ix[
                    i, "Noise_Gain_Lin"] + self.stack.ix[i, "Noise_added_Lin"]
                self.stack.ix[i, "Noise_In_Lin_Hz"] = self.stack.ix[i - 1, "Noise_Out_Lin_Hz"]

        self.stack["Eq_Noise_Out_dBm_Hz"] = 10*np.log10(self.stack["Eq_Noise_Out_Lin_Hz"]/.001)
        self.stack["Eq_Noise_In_dBm_Hz"] = 10*np.log10(self.stack["Eq_Noise_In_Lin_Hz"]/.001)

        self.stack["Noise_Out_dBm_Hz"] = 10 * np.log10(self.stack["Noise_Out_Lin_Hz"] / .001)
        self.stack["Noise_In_dBm_Hz"] = 10 * np.log10(self.stack["Noise_In_Lin_Hz"] / .001)

    def populate_vswr(self,vswr):

        space = len(self.stack['VSWR'])
        values = len(vswr)

        if space == values:
            self.stack['VSWR'] = vswr
        else:
            raise ValueError('Error(populate VSWR): there are ', space, ' spaces in the array and you gave me ', values, ' values for VSWR')

        if bin(self.flags & self.VSWR_mask).count("1") == 0:
            self.flags = self.flags + self.VSWR_mask


    def calculate_art(self):

        temp = self.stack.copy()
        temp_gain = temp['Gain_dB']

        tau = np.zeros_like(temp_gain)
        self.stack['|aRT|_lin'] = np.zeros_like(temp_gain)
        self.stack['|acbl|max_lin'] = np.zeros_like(temp_gain)
        self.stack['|acbl|min_lin'] = np.zeros_like(temp_gain)
        self.stack['|acbl|max_dB'] = np.zeros_like(temp_gain)
        self.stack['|acbl|min_dB'] = np.zeros_like(temp_gain)

        for i in range(0, len(self.stack["Type"])):

            if self.stack.ix[i, "Type"] == 'cable' or self.stack.ix[i, "Type"] == 'tx_line':
                tau = 10**(self.stack.ix[i, 'Gain_dB']/20)
                self.stack.ix[i,'|aRT|_lin'] = (tau**2)*((self.stack.ix[i-1, 'VSWR']-1)/(self.stack.ix[i-1, 'VSWR']+1))*((self.stack.ix[i, 'VSWR']-1)/(self.stack.ix[i, 'VSWR']+1))
                self.stack.ix[i,'|acbl|max_lin'] = tau  / (1 - self.stack.ix[i,'|aRT|_lin'])
                self.stack.ix[i,'|acbl|min_lin'] = tau /  (1 + self.stack.ix[i, '|aRT|_lin'])

                self.stack.ix[i, '|acbl|max_dB'] = 20*np.log10(self.stack.ix[i,'|acbl|max_lin'])
                self.stack.ix[i, '|acbl|min_dB'] = 20*np.log10(self.stack.ix[i,'|acbl|min_lin'])


        print(tau)



if __name__ == '__main__':

    gains = [    20,      -1.25,      -0.5,           -2,     -10,   10,     0]#np.random.uniform(-15,26,4)
    NF    = [     3,       1.25,       0.5,            2,      10,    5,    30]#np.random.uniform(0.5,26,4)
    ITOI  = [    20,        100,       100,          100,      16,   10,    24]#[8,11.2,np.nan,1]#np.random.uniform(-10,26,4)
    IP2   = [    30,        100,       100,          100,      26,   20,    34]
    inputs= [     1,          1,         1,            1,       1,    1,     1]
    IP1dB  =[    10,        100,       100,          100,       6,    0,     0]
    VSWR =  [   1.8,        2.0,       2.1,          2.0,     3.1,  6.1,   1.1]
    Stage_Bandwidth = 50e6*np.ones_like(gains)
    names = ['amp',     'cable',  'filter',     'cable', 'mixer_wb','vga', 'ADC']#,'split','mixer', 'filter','lin','filter']#['filter','thing2','thing3','thing4','thing5','thing7','mixer']


    rf_stackup = RF_stackup(names)
    rf_stackup.stack['IP1dB'] = IP1dB
    rf_stackup.stack['Device Inputs'] = inputs
    rf_stackup.populate_gains_dB(gains)
    rf_stackup.populate_vswr(VSWR)
    rf_stackup.calculate_art()
    rf_stackup.calculate_OP1dB_dB()
    rf_stackup.cascade_OP1dB_lin()
    rf_stackup.calculate_Coherent_Gain(inputs)
    rf_stackup.calculate_siginout(-20)
    rf_stackup.populate_NF_dB(NF)
    rf_stackup.cascade_NF()
    rf_stackup.calculate_noise_gain()

    rf_stackup.populate_ITOI_dB(ITOI)
    rf_stackup.calculate_OTOI_dB()
    #rf_stackup.calculate_lin_OTOI()
    rf_stackup.cascade_OTOI()

    rf_stackup.populate_IP2_dB(IP2)
    rf_stackup.calculate_OP2_dB()
    rf_stackup.cascade_OP2()


    print(rf_stackup.flags)
    print(rf_stackup.stack)
    print(np.shape(rf_stackup.stack))


    #then noise input noise power, input power, elements or combiner networks (for summing in arrays)
    #then vswr inputs, then fillnans on vswr, max_gain _min_gain, sigma,

