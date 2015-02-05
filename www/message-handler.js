// This recieves messages of type "restore.like.state" from the server.
Shiny.addCustomMessageHandler("restore.like.state",
	function(message) {
		for (property in message) {
			var obj = document.getElementById(property);
			if (obj != null) {
				console.log(property)
				//console.log(obj)
			}
		}
		//console.log(message);
		//console.log(JSON.stringify(message));
	}
);
