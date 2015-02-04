SELECT `instance.pid`, `instance.time`, `session.id`, `session.counter`,
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
