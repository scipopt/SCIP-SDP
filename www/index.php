<!DOCTYPE html>
<html>
<head>
   <link href='http://fonts.googleapis.com/css?family=Open+Sans' rel='stylesheet' type='text/css'>
   <meta http-equiv="content-type" content="text/html; charset=UTF-8">
   <meta name="viewport" content="width=device-width, initial-scale=1.0, maximum-scale=1, user-scalable=no">
   <title>SCIP-SDP</title>
   <!--   load css styles -->
   <link rel="stylesheet" type="text/css" href="bootstrap/css/bootstrap.min.css">
   <link rel="stylesheet" type="text/css" href="bootstrap/css/custom.css">
   <link rel="stylesheet" type="text/css" href="bootstrap/css/custom-sdp.css">
   <!-- jQuery (necessary for Bootstrap's JavaScript plugins) -->
   <script src="bootstrap/js/jquery.js"></script>
   <!-- Include all compiled plugins (below), or include individual files as needed -->
   <script src="bootstrap/js/custom.js"></script>
   <script src="bootstrap/js/bootstrap.min.js"></script>
</head>
<body>
   <?php include('header.inc'); ?>
   <?php include('sdp-banner.inc'); ?>
   <div class="container bs-docs-container">
      <div class="row">
         <div class="col-md-3">
            <?php include('sdp-menu.inc'); ?>
         </div>
         <!-- end col-md-3 -->
         <div class="col-md-9">
            <div>
               <h1 id="about">About</h1>
               <?php include('about.inc'); ?>
            </div>
            <div>
               <h1 id="news">News</h1>
               <?php include('news.inc'); ?>
            </div>
            <div>
               <h1 id="license">License</h1>
               <?php include('license.inc'); ?>
            </div>
            <div>
               <h1 id="bugs">Bugs</h1>
               <?php include('bugs.inc'); ?>
            </div>
            <div>
               <h1 id="developers">Developers</h1>
               <?php include('developers.inc'); ?>
            </div>
            <div>
               <h1 id="download">Download</h1>
               <?php include('download.inc'); ?>
            </div>
         </div>
         <!-- end col-md-9 -->
      </div>
      <!-- end row -->
   </div>
   <!-- end container -->
   <footer class="bs-footer" role="contentinfo">
      <div class="container">
         <div class="row">
            <div class="col-md-12">
                <?php include('footer.inc'); ?> 
            </div>
         </div>
      </div>
   </footer>
</body>
</html>

