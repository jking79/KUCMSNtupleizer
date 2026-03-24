#change indices for what the substrs should be based on whatever the input name is
def uniformize_sig_name(instr):
    path = instr.split("/")[:-1]
    path = "/".join(path)
    splitstr = instr[len(path)+1:].split("_")
    evtfilter = splitstr[2]
    ver = splitstr[4]
    sms_proc = splitstr[1]
    data_type = splitstr[6]
    mgl = splitstr[-5]
    mn2 = splitstr[-4]
    mn1 = splitstr[-3]
    ct = splitstr[-2]
    if "p" not in ct:
        ct = "ct0p"+ct[2:]
    newname = f"{path}/SMS_{evtfilter}_{ver}_{sms_proc}_{data_type}_{mgl}_{mn2}_{mn1}_{ct}_rjrskim.root"
    return newname
