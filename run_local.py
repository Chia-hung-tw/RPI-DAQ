import datetime,ctypes,yaml
import rpi_daq, unpacker
import skiroc2cms_bit_string as sk2conf
from optparse import OptionParser

class yaml_config:
    yaml_opt=yaml.YAMLObject()
    
    def __init__(self,fname="default-config.yaml"):
        with open(fname) as fin:
            self.yaml_opt=yaml.safe_load(fin)
            
    def dump(self):
        return yaml.dump(self.yaml_opt)

    def dumpToYaml(self,fname="config.yaml"):
        with open(fname,'w') as fout:
            yaml.dump(self.yaml_opt,fout)
        
def get_comma_separated_args(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

if __name__ == "__main__":
    
    parser = OptionParser()
    parser.add_option("-a", "--compressRawData",dest="compressRawData",action="store_true",
                      help="option to compress (ie remove the '1000' before each word) raw data",default=False)
    
    parser.add_option("-b", "--moduleNumber", dest="moduleNumber",type="int",action="store",
                      help="moduleNumber", default=63)
    
    parser.add_option("-c", "--nEvent", dest="nEvent",type="int",
                      help="number of event",default=100)

    parser.add_option("-d", "--externalChargeInjection", dest="externalChargeInjection",action="store_true",
                      help="set to use external injection",default=False)
    parser.add_option("-e", "--acquisitionType", dest="acquisitionType",choices=["standard","sweep","fixed","const_inj"],
                      help="method for injection", default="standard")
    parser.add_option('-f', '--channelIds', dest="channelIds",action="callback",type=str,
                      help="channel Ids for charge injection", callback=get_comma_separated_args, default=[])
    parser.add_option('-g','--injectionDAC',dest="injectionDAC",type="int",action="store",default=1000,
                      help="DAC setting for injection when acquisitionType is const_inj")
    (options, args) = parser.parse_args()
    print(options)

    conf=yaml_config()
    conf.yaml_opt['daq_options']['compressRawData']=options.compressRawData
    conf.yaml_opt['daq_options']['nEvent']=options.nEvent
    conf.yaml_opt['daq_options']['acquisitionType']=options.acquisitionType
    conf.yaml_opt['daq_options']['externalChargeInjection']=options.externalChargeInjection
    conf.yaml_opt['daq_options']['injectionDAC']=options.injectionDAC
    for i in options.channelIds:
        conf.yaml_opt['daq_options']['channelIds'].append(int(i))

    conf.yaml_opt['glb_options']['moduleNumber'] = options.moduleNumber
    daq_options=conf.yaml_opt['daq_options']
    glb_options=conf.yaml_opt['glb_options']

    
    print "daq options = "+yaml.dump(daq_options)
    print "Global options = "+yaml.dump(glb_options)

    the_bit_string=sk2conf.bit_string()
    the_bit_string.Print()
    if daq_options['externalChargeInjection']==True:
        if len(daq_options['channelIds'])>0:
            the_bit_string.set_channels_for_charge_injection(daq_options['channelIds'])
        else:
            print("Option channelIds should not be empty if charge injection is set")        

    if len(daq_options['channelIdsToMask'])>0:
        the_bit_string.set_channels_to_mask(daq_options['channelIdsToMask'])
        
    if len(daq_options['channelIdsDisableTOT'])>0:
        the_bit_string.set_channels_to_disable_trigger_tot(daq_options['channelIdsDisableTOT'])

    if len(daq_options['channelIdsDisableTOA'])>0:
        the_bit_string.set_channels_to_disable_trigger_toa(daq_options['channelIdsDisableTOA'])


    the_bit_string.Print()
    the_bits_c_uchar_p=the_bit_string.get_48_unsigned_char_p()
    print( [hex(the_bits_c_uchar_p[i]) for i in range(48)] )


    the_time=datetime.datetime.now()
    if glb_options['storeYamlFile']==True:
        yamlFileName=glb_options['outputYamlPath']+"/Module"+str(glb_options['moduleNumber'])+"_"
        yamlFileName=yamlFileName+str(the_time.day)+"-"+str(the_time.month)+"-"+str(the_time.year)+"_"+str(the_time.hour)+"-"+str(the_time.minute)
        yamlFileName=yamlFileName+".yaml"
        print("\t save yaml file : ",yamlFileName)
        conf.dumpToYaml(yamlFileName)
    
    rawFileName=glb_options['outputRawDataPath']+"/Module"+str(glb_options['moduleNumber'])+"_"
    rawFileName=rawFileName+str(the_time.day)+"-"+str(the_time.month)+"-"+str(the_time.year)+"_"+str(the_time.hour)+"-"+str(the_time.minute)
    rawFileName=rawFileName+".raw"
    print("\t open output file : ",rawFileName)
    logfile = open("runinfo.log",'w')
    logfile.write(rawFileName)
    logfile.close()

    outputFile = open(rawFileName,'wb')

    theDaq=rpi_daq.rpi_daq(daq_options)
    outputBitString=theDaq.configure(the_bits_c_uchar_p)
    print("\t write bits string in output file")
    byteArray = bytearray(outputBitString)
    outputFile.write(byteArray)

    data_unpacker=unpacker.unpacker(daq_options['compressRawData'])
    #start_dac = raw_input("Start of the dac value: ")
    #spacing_dac = raw_input("Spacing of the dac value: ")
    
    for event in range(daq_options['nEvent']):
        #dac = int(start_dac) + int(spacing_dac) * event
        dac = 0
        if dac > 4095:
            dac = 4095
        rawdata=theDaq.processEvent(dac)
        #data_unpacker.unpack(rawdata)
        #data_unpacker.showData(event)
    
        byteArray = bytearray(rawdata)
        outputFile.write(byteArray)
        
        if event%10==0:
            print("event number ",event)