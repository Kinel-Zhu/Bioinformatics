##方法3个： sqldf RMySQL RODBC
# library(sqldf)
# sqldf("select * from drug",dbname="test",drv="MySQL",user="root",password="Z123",host="127.0.0.1",port=3306)
# sqldf("select * from t_data limit 0,10")

##install.packages("RMySQL")
##library(RMySQL)
##conn <- dbConnect(MySQL(), dbname = "cancer", username="root", password="Z123", host="127.0.0.1", port=3306)
##clin<-dbReadTable(conn, "cancer_sample_inf") 

###install.packages("RODBC")
###library(RODBC) 使用与sqlsever

##选择RMySQL
library(RMySQL)
conn <- dbConnect(MySQL(), dbname = "cancer", username="root", password="Z123", host="127.0.0.1", port=3306)

dbListTables(conn)##查表
clin<-dbReadTable(conn, "cancer_sample_inf")##读表

res1 <- dbSendQuery(conn, "SELECT  cancer_name FROM cancer_sample_inf")##语句查询，没有数据取出来
data <- dbFetch(res1, n=2)##取前两条数据 ，无n属性则取全部数据

res2<-dbGetQuery(conn, "SELECT death_num FROM cancer_sample_inf where cancer_name='BLCA' ")##取数据的第二种方法,并条件筛选
res3<-dbGetQuery(conn, "SELECT death_num FROM cancer_sample_inf")

dbWriteTable(conn, "tablename", data)##写数据
dbDisconnect(conn) ##关闭连接
