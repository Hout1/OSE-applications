#!/usr/bin/env python
# coding: utf-8

# In[12]:


#!/usr/bin/env python
# coding: utf-8

# In[1]:

def main():

    from datetime import datetime, timedelta, date
    from configparser import ConfigParser
    import os
    import schedule
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart #python3
    from email.mime.base import MIMEBase #python3
    from email import encoders #python3

    def send_mail():
        msg = MIMEMultipart()
        msg['From'] = 'wfp.hq.dbtrack@gmail.com'
        msg['Subject'] = "Flooding Event Alert - BETA version"

        body = """<div style="color:rgb(28, 130, 196)"><font size="+2"><b><img id="myimage" src="https://mw1.google.com/crisisresponse/icons/un-ocha/disaster_flood_100px_icon.png"></img>
        <br />Flooding Alert!</b> </font></div>
        
        <br/><b>Data Availability for Dashboard Generation</b></b>
        <br />
        <br/>MODIS (Global) Automatically generated</b>
        <br/>AER (Africa) Manually launched by HQGIS:</b>
        <br/>Sentinel (Global) Manually generated</b>
        <br/>Flood impact assessment: Manually launched by HQGIS</b>
        
        <br/>Please see attached flood summary and initial impact assessment.</b>
        <br />
        <br/><b>Emergency Mapping/Charter Activation Links:</b></b>
        <br />
        <br/><b>--------------------------------------------------------------------------------</b></b>
        <br/>This is an automatically generated email, please do not reply.<br/>
        <br/>Service provided by 
        <br/><i>UNITED NATIONS WORLD FOOD PROGRAMME
        <br/>WFP Emergency Preparedness & Support Response Division
        <br/>Contact: <a href="mailto:hq.gis@wfp.org">HQ Geospatial Support Unit</a></i></p>
        """
        
        alert_fn = alert_pdf_path + alert_data_name

        alerts = MIMEBase('application', "octet-stream")

        alerts.set_payload(open(alert_fn, "rb").read())

        encoders.encode_base64(alerts)

        alerts.add_header('Content-Disposition', 'attachment; filename=%s'%alert_data_name)

        msg.attach(alerts)

        msg.attach(MIMEText(body, 'html'))

    #     s = smtplib.SMTP('smtp.gmail.com', 587) #Anywhere
    #     s.starttls() #Anywhere

        s = smtplib.SMTP_SSL('smtp.gmail.com', 465) # WFP domain

        s.ehlo()
        s.login(smtp_un, smtp_pw)
        sender = 'wfp.hq.dbtrack@gmail.com'

        recipients = ['michaelandrew.manalili@gmail.com']#,'michael.manalili@wfp.org','abdel-lathif.younous@wfp.org','rohini.swaminathan@wfp.org']
        #recipients = ['rohini.swaminathan@wfp.org','sirio.modugno@wfp.org','stefano.cairo@wfp.org','michael.manalili@wfp.org','abdel-lathif.younous@wfp.org','michaelandrew.manalili@gmail.com']    

        s.sendmail(sender, recipients, str(msg))
        s.quit()
        print('FLOOD ALERT sent to subscribers!')

    now = datetime.today()
    d = str(now)
    date = d[0:10]
    date_cut = date.replace('-', '')

    script = os.path.dirname(os.path.realpath("__file__"))
    root = os.path.dirname(script)
    
    alert_pdf_path = root + '/alert/'
    
    alert_data_name = 'WFP_FloodSummary_' + date_cut + '_merged' + '.pdf'

    config = ConfigParser()
    config.read(os.path.join(root, "config.txt"))

    config.read(os.path.join(root, "config.txt"))
    smtp_un = config.get('dbtrack','dbt_uname')
    smtp_pw = config.get('dbtrack','dbt_pword')
    print('config loaded.')

    send_mail()
    
if __name__ == "__main__":
    main()


# In[ ]:




