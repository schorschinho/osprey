function[] = OspreyAutoMail(MRSCont)
% function[MRSCont] = OspreyAutoMail(MRSCont)
% Function for emailing report HTMLs to a recipient list.
% 
% Input:    MRSCont=Osprey container with appropriate options for AutoMail.
%
% Requires specification of:
%   1) A JSON, MRSCont.opts.mailto.config, containing source email and
%      associated password.
%   2) Recipient list, MRSCont.opts.mailto.recipients, as a cell array.
% 
% Currently supports:
%   1) Non-authenticated email servers.
%   2) Gmail adresses configured to allow MATLAB as external app.
%
% C.W.Davies-Jenkins, Johns Hopkins Medicine 2022

% Read config JSON:
fid = fopen(MRSCont.opts.mailto.config); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
if strcmp('win',osp_platform('filesys'))
    str = strrep(str,'\','\\');
end
str = replace(str, whitespacePattern + '"', '"');
Config  = jsondecode(str);

% Set source email and password using config json
setpref('Internet','E_mail',Config.SourceEmail);
setpref('Internet','SMTP_Username',Config.SourceEmail);
setpref('Internet','SMTP_Password',Config.Password);

% % Set preferences for specific email server
Src = lower(extractAfter(Config.SourceEmail,'@'));
switch Src
    case 'gmail.com'
        %set up SMTP service for Gmail
        setpref('Internet','SMTP_Server','smtp.gmail.com');
        % Gmail server.
        props = java.lang.System.getProperties;
        props.setProperty('mail.smtp.auth','true');
        props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
        props.setProperty('mail.smtp.socketFactory.port','465');
    otherwise
        warning('Email source: "%s" is not currently explicitly supported. Trying anyway!...',Src)
end

% Look for relevant files to include in email
Include = cell(0);
if exist(fullfile(MRSCont.outputFolder,'Reports'),'dir') % If reports exist, send them
    Include = [Include,fullfile(MRSCont.outputFolder,'Reports/*.html')];
    NReports = length(dir(fullfile(MRSCont.outputFolder,'Reports/*.html')));
else
    NReports = 0;
end
if exist(fullfile(MRSCont.outputFolder, 'LogFile.txt'),'file') % If log file exists, send it
    Include = [Include,fullfile(MRSCont.outputFolder, 'LogFile.txt')];
end
zip(fullfile(MRSCont.outputFolder,'MRS_Reports'), Include); % Create zipfile with includes

Subject = 'Osprey auto-email'; % Email subject
Message = sprintf('Hi,\n\nHere are your MRS results!\n We found %i reports to attach.\n\nBest wishes,\nTeam Osprey.',NReports); % main body of email.

sendmail(MRSCont.opts.mailto.recipients, Subject, Message, fullfile(MRSCont.outputFolder,'MRS_Reports.zip')); % Send email with zipfile

% Cleanup:
delete(fullfile(MRSCont.outputFolder,'MRS_Reports.zip')); % Remove zipfile

% Remove the preferences (for privacy reasons)
setpref('Internet','E_mail','');
setpref('Internet','SMTP_Server','''');
setpref('Internet','SMTP_Username','');
setpref('Internet','SMTP_Password','');

end