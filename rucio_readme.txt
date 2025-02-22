Rucio
Creating a rule

To make a basic rule to copy a CMS dataset or block to a specific site, the easiest method is to use the Rucio Web UI via the “Request New Rule” option under the “Data Transfers” tab
￼
or using this link: (https://cms-rucio-webui.cern.ch/r2d2/request
￼
). The UI will have you fill out each field in order to create a rule, namely:
* The DID for which you are trying to create a rule (either a container or dataset). The scope must precede the name of the DID (usually this will be “cms:”). The Web UI includes the option to create a rule for several DIDs at once; to do this, choose the “List of DIDs” option before entering a DID name.
* Choosing the RSE the rule shall be associated with. The system will check the storage available to your account at the specified RSE. If you do not have any available quota, you can choose the "I want to ask for approval" option to create the rule and alert the site admin.
* All extra options, e.g., lifetime, number of copies, and whether the rule shall be created asynchronously.
To create a rule using the CLI, follow the template below: rucio add-rule cms:/CMS/DATA/SET/NAME 1 T2_MY_SITE or rucio add-rule cms:/CMS/DATA/SET/NAME#BLOCK-NAME 1 T2_MY_SITE Keep in mind the 'cms:' is required, this is what Rucio calls a 'scope'. All centrally produced CMS data is in the 'cms' scope. Also note that what CMS calls a block, Rucio calls a 'DATASET' and what CMS calls a dataset, Rucio calls a 'CONTAINER'. If dataset is huge please add --asynchronous it will allows create rule asynchronously. Rucio rules are very flexible, so rucio add-rule --help will give you more information on what can be done. This command returns a unique rule ID which you can save and reference later to query the status of a rule.
However, before making a rule, you should know if you have sufficient quota for the desired RSE (or RSE expression). You can find out by following RucioUserDocsQuotas. If you do not have a quota, or you have reached your quota limit, you can exceptionally make a rule request which notifies the site administrator of your request and allows them to approve the rule irrespective of your quota. To do this, add the option '--ask-approval' to the 'add-rule' command above, e.g. rucio add-rule --ask-approval cms:/CMS/DATA/SET/NAME 1 T2_MY_SITE Your request is more likely to be approved if you also specify a --lifetime (in seconds; for reference, 30 days is 2592000 seconds) that you require the data to be protected by the rule, as this lets site administrators know how long the data will occupy space at their site.
If you have a rule holding data at a given location, it will not disappear from that location. If you don't have a rule, but data presently exists at some location, there is no guarantee it will stay. Each rule you make contributes to your account's RSE usage, regardless of whether or not the data is protected by other rules, so you should be sure to request enough quota to hold all your rules or make each rule with --ask-approval and a --lifetime.

Querying rules
* rucio list-rules --account $RUCIO_ACCOUNT List all the rules you have and their state
* rucio rule-info [RULE_HASH] monitor the progress of your rule and any transfers it may have initiated
* Using the Web UI, a rule can be queried directly using rule ID, returning info equivalent to the rule-info function. To check your own rules, navigate to the “List my rules” section under the “Data Transfer” tab.

rucio list-rules --account jaking
rucio delete-rule #
rucio rule-info #


Removing rules #
If you created a rule with an explicit lifetime, it will automatically get cleaned up after it expires. If you created a rule without a lifetime, you may explicitly delete it via the command c [RULE_HASH]. When you delete the rule, your quota will be returned. The data replica may or may not be deleted from the site, depending on the actual storage usage at the site and any other users' rules that apply to the same data at the same site.
rucio list-rules --account jaking
if you intend to use a Rucio rule to “subscribe” files to your site, you can do this for just one DBS block or even a single file.
See e.g. rucio add-rule --help A DID could be e.g. cms:<LogicalFileName> or cms:<DBS block name> If you want/need a different mix, you need to create a container in your user scope inclluding CMS did’s as proper.
Or you can catch files “on the fly” after CRAB recalled them for you (they should linger on disk for months) and simply copy them, e.g. via rucio get or via gfal-copy if you prefere, in the latter case you may want to use rucio list-file-replicas cms:<LFN> to get the PFN to use as source.
