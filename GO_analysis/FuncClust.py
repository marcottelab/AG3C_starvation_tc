#!python
# by courtesy of HuangYi @ 20110424

def DAVIDenrich(listF, idType, bgF='', resF='', bgName = 'Background1',listName='List1', category = '', thd=0.1, ct=2):
        from suds.client import Client
        import os
	try:
	   
	    overlap=3
	    initalSeed=3 
	    finalSeed=3
	    linkage=0.5
	    kappa=1
	    if len(listF) > 0 and os.path.exists(listF):
		inputListIds = ','.join(open(listF).read().split('\n'))
		print 'List loaded.'        
	    else:
		print 'No list loaded.'
		raise

	    flagBg = False
	    if len(bgF) > 0 and os.path.exists(bgF):
		inputBgIds = ','.join(open(bgF).read().split('\n'))
		flagBg = True
		print 'Use file background.'
	    else:
		print 'Use default background.'

	    client = Client('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl',faults=True)
	    print 'User Authentication:',client.service.authenticate('johnrhouser@utexas.edu')

	    listType = 0
	    print 'Percentage mapped(list):', client.service.addList(inputListIds,idType,listName,listType)
	    if flagBg:
		listType = 1
		print 'Percentage mapped(background):', client.service.addList(inputBgIds,idType,bgName,listType)

	    print 'Use categories:', client.service.setCategories(category)
	    #chartReport = client.service.getChartReport(thd,ct)
	    chartReport=client.service.getTermClusterReport(overlap,initalSeed,finalSeed,linkage,kappa)
	    #chartReport=client.service.getGeneClusterReport(overlap,initalSeed,finalSeed,linkage,kappa)

	    chartRow = len(chartReport)
	    print 'Total chart records:',chartRow
	    
	    if len(resF) == 0 or not os.path.exists(resF):
		if flagBg:
		    resF = listF + '.withBG.chartReport'
		else:
		    resF = listF + '.clustReport'
	    with open(resF, 'w') as fOut:
		fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
		for row in chartReport:
		    rowDict = dict(row)
		    for j in range(len(rowDict['simpleChartRecords'])):
			    categoryName = str(rowDict['simpleChartRecords'][j]['categoryName'])
			    termName = str(rowDict['simpleChartRecords'][j]['termName'])
			    listHits = str(rowDict['simpleChartRecords'][j]['listHits'])
			    percent = str(rowDict['simpleChartRecords'][j]['percent'])
			    ease = str(rowDict['simpleChartRecords'][j]['ease'])
			    Genes = str(rowDict['simpleChartRecords'][j]['geneIds'])
			    listTotals = str(rowDict['simpleChartRecords'][j]['listTotals'])
			    popHits = str(rowDict['simpleChartRecords'][j]['popHits'])
			    popTotals = str(rowDict['simpleChartRecords'][j]['popTotals'])
			    foldEnrichment = str(rowDict['simpleChartRecords'][j]['foldEnrichment'])
			    bonferroni = str(rowDict['simpleChartRecords'][j]['bonferroni'])
			    benjamini = str(rowDict['simpleChartRecords'][j]['benjamini'])
			    FDR = str(rowDict['simpleChartRecords'][j]['afdr'])
			    rowList = [categoryName,termName,listHits,percent,ease,Genes,listTotals,popHits,popTotals,foldEnrichment,bonferroni,benjamini,FDR]
			    fOut.write('\t'.join(rowList)+'\n')
		
		    fOut.write('\n')
		print 'write file:', resF, 'finished!'

	except:
		print "Warning Something went wrong in functional clustering of" + str(listF) + "perhaps there aren't enough genes in this list"
if __name__ == '__main__':
    DAVIDenrich(listF = './RNA_analysis/RNA_all_combined_conv_list.txt', idType = 'ENTREZ_GENE_ID',bgName='Escherichia coli', category = 'GOTERM_BP_FAT')
    #DAVIDenrich(listF = './list1.txt', idType = 'AFFYMETRIX_3PRIME_IVT_ID', listName = 'list1', category = 'GOTERM_BP_FAT')
        
