# Bloomberg CDS Data

The CDS data are fetched using Bloomberg Excel plug-ins. The security string is of the structure `[$COUNTRY_NAME] CDS USD [$TENOR]Y D14`.
A few note-worthy things:
1. Many countries not in original AEJ list have little to no data, especially on 1Y CDS. e.g. UK has 1 data point, which is highly implausible. LPPS also explained that their choice of countries are limited by the availability of Bloomberg data. 

# Discount Value

1. download Treasury constant maturity data [here](https://www.federalreserve.gov/datadownload/Choose.aspx?rel=H15)