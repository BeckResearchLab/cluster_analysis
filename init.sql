DROP TABLE IF EXISTS log;
CREATE TABLE log (
	id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
	`instance.pid`	INT,
	`instance.time`	BIGINT(20),
	`session.id`		INT,
	`session.counter` INT,
	k				INT,
	kLikeReason		VARCHAR(128),
	cluster			INT,
	clusterLikeReason	VARCHAR(128),
	clusterProfilePlotSampleNames VARCHAR(32),
	clusterProfilePlotTracks		VARCHAR(32),
	clusterDisplayMotif1GeneProfile	INT,
	clusterDisplayMotif2GeneProfile	INT,
	clusterDisplayMotif4GeneProfile	INT,
	clusterDisplayMotif3GeneProfile	INT,
	clusterMotif1LikeReason	VARCHAR(128),
	clusterMotif2LikeReason	VARCHAR(128),
	clusterMotif3LikeReason	VARCHAR(128),
	clusterMotif4LikeReason	VARCHAR(128),
	clusterSelectedRows	TEXT,
	myClusterGenes		TEXT,
	myClusterRecruitN	INT,
	myClusterRecruitBy	VARCHAR(32),
	myClusterLikeReason	VARCHAR(128),
	myClusterProfilePlotTracks	TEXT,
	myClusterProfilePlotSampleNames	TEXT,
	myClusterDisplayMotif1GeneProfile INT,
	myClusterDisplayMotif2GeneProfile INT,
	myClusterDisplayMotif3GeneProfile INT,
	myClusterDisplayMotif4GeneProfile INT,
	myClusterMotif1LikeReason	VARCHAR(128),
	myClusterMotif2LikeReason	VARCHAR(128),
	myClusterMotif3LikeReason	VARCHAR(128),
	myClusterMotif4LikeReason	VARCHAR(128),
	myClusterSelectedRows	TEXT,
	searchText		TEXT,
	clusterSearchResultSelectedRow	TEXT,
	`blastn.input`	TEXT,
	blastnDatabase	VARCHAR(32),
	`blastp.input`	TEXT,
	row_names		CHAR(1),	# this is a field added by R, we ignore it
	INDEX(kLikeReason),
	INDEX(clusterLikeReason),
	INDEX(clusterMotif1LikeReason),
	INDEX(clusterMotif2LikeReason),
	INDEX(clusterMotif3LikeReason),
	INDEX(clusterMotif4LikeReason),
	INDEX(myClusterLikeReason),
	INDEX(myClusterMotif1LikeReason),
	INDEX(myClusterMotif2LikeReason),
	INDEX(myClusterMotif3LikeReason),
	INDEX(myClusterMotif4LikeReason)
);


DROP VIEW IF EXISTS log_likes;
CREATE VIEW log_likes AS
SELECT id,
			CASE
				WHEN kLikeReason <> "" THEN "k"
				WHEN clusterLikeReason <> "" THEN "cluster"
				WHEN clusterMotif1LikeReason <> "" THEN "clusterMotif1"
				WHEN clusterMotif2LikeReason <> "" THEN "clusterMotif2"
				WHEN clusterMotif3LikeReason <> "" THEN "clusterMotif3"
				WHEN clusterMotif4LikeReason <> "" THEN "clusterMotif4"
				WHEN myClusterLikeReason <> "" THEN "myCluster"
				WHEN myClusterMotif1LikeReason <> "" THEN "myClusterMotif1"
				WHEN myClusterMotif2LikeReason <> "" THEN "myClusterMotif2"
				WHEN myClusterMotif3LikeReason <> "" THEN "myClusterMotif3"
				WHEN myClusterMotif4LikeReason <> "" THEN "myClusterMotif4"
				ELSE NULL
			END
			AS liked_how, 
			CASE
				WHEN kLikeReason <> "" THEN kLikeReason
				WHEN clusterLikeReason <> "" THEN clusterLikeReason
				WHEN clusterMotif1LikeReason <> "" THEN clusterMotif1LikeReason
				WHEN clusterMotif2LikeReason <> "" THEN clusterMotif2LikeReason
				WHEN clusterMotif3LikeReason <> "" THEN clusterMotif3LikeReason
				WHEN clusterMotif4LikeReason <> "" THEN clusterMotif4LikeReason
				WHEN myClusterLikeReason <> "" THEN myClusterLikeReason
				WHEN myClusterMotif1LikeReason <> "" THEN myClusterMotif1LikeReason
				WHEN myClusterMotif2LikeReason <> "" THEN myClusterMotif2LikeReason
				WHEN myClusterMotif3LikeReason <> "" THEN myClusterMotif3LikeReason
				WHEN myClusterMotif4LikeReason <> "" THEN myClusterMotif4LikeReason
			END
			AS liked_reason
	FROM log
	WHERE 
			CASE
				WHEN kLikeReason <> "" THEN 1
				WHEN clusterLikeReason <> "" THEN 1
				WHEN clusterMotif1LikeReason <> "" THEN 1
				WHEN clusterMotif2LikeReason <> "" THEN 1
				WHEN clusterMotif3LikeReason <> "" THEN 1
				WHEN clusterMotif4LikeReason <> "" THEN 1
				WHEN myClusterLikeReason <> "" THEN 1
				WHEN myClusterMotif1LikeReason <> "" THEN 1
				WHEN myClusterMotif2LikeReason <> "" THEN 1
				WHEN myClusterMotif3LikeReason <> "" THEN 1
				WHEN myClusterMotif4LikeReason <> "" THEN 1
				ELSE 0
			END
			= 1
	ORDER BY `instance.time`, `instance.pid`, `session.id`, `session.counter`
;
