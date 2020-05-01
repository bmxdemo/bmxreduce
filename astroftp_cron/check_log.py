receivers = [line.strip() for line in 
            open('/astro/u/bmx/bmxreduce/astroftp_cron/emailtos').readlines()]

import datetime
import smtplib, ssl


def check_file(fname):
    lines = open(fname).readlines()[:5]
    line0, line1 = lines[:2]
    exptime = datetime.datetime.now().strftime("%Y-%m-%d %H.05")
    lasttime = line0[:len(exptime)]
    ok = False
    if exptime != lasttime:
        msg = ("Last time record seems wrong (dead computer?) \n      Expected: %s\n      Got: %s"%
                 (exptime, lasttime))
    else:
        def getstatus(line):
            return line.split(":")[1].split()[0]
        status0,status1 = [getstatus(line) for line in [line0,line1]]
        if getstatus(line0) != getstatus(line1):
            msg = "Status change from %s to %s"%(status1,status0)
        else:
            msg = "No change"
            ok = True

    return ok, msg, lines




ok1, msg1, lines1 = check_file("/gpfs02/astro/www/bmx/daqlog/bmx1.log")
ok2, msg2, lines2 = check_file("/gpfs02/astro/www/bmx/daqlog/bmx2.log")
if (not ok1) or (not ok2):
    print ("Sending emails...")
    port = 465  # For SSL
    smtp_server = "rcf.rhic.bnl.gov"
    sender_email = "anze@bnl.gov"  # Enter your address
    password = open('/astro/u/bmx/bmxreduce/astroftp_cron/emailpw').readlines()[0].replace('\n','')
    message = """\
From: BMX
Subject: BMX alert

This is an automated email from a script living at
/astro/u/bmx/bmxreduce/astroftp_cron/check_log.py 

BMX1: %s
Last 5 entries:
%s
BMX2: %s
Last 5 entries:
%s"""%(msg1, "".join(lines1), msg2, "".join(lines2))

    context = ssl.create_default_context()
    with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
        server.login(sender_email, password)
        for receiver_email in receivers:
            server.sendmail(sender_email, receiver_email, message)
else:
    print ("Everything OK...")
