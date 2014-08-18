/*
VCell-ID is is a graphical user interface which provides easy testing options
that help to choose the correct parameters to process microscope fluorescence
images with CellID (a software developped by Andrew Gordon to identify cells in 
images and to calculate a number of statistics, including statistics derived from corresponding fluorescence images).
Copyright (C) 2008 Ariel Chernomoretz

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
***************************************************************************
/* VCellID makes use of gtkimageviewer library written by
 * Björn Lindqvist <bjourne@gmail.com>, and distributed under
 * the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2, or
 * (at your option) any later version.

**********************************************************************
Start-copyright-notice-for-libtiff
Libtiff software is used by Cell-ID for some of the reading in of
TIF image file data and also for creating new TIF files. Libtiff is
available at http://www.remotesensing.org/libtiff/. The libtiff software
was written by Sam Leffler while working for Silicon Graphics and
contains the following copyright notice:

   "Copyright (c) 1988-1997 Sam Leffler
    Copyright (c) 1991-1997 Silicon Graphics, Inc.

    Permission to use, copy, modify, distribute, and sell this software and 
    its documentation for any purpose is hereby granted without fee, provided
    that (i) the above copyright notices and this permission notice appear in
    all copies of the software and related documentation, and (ii) the names
    of Sam Leffler and Silicon Graphics may not be used in any advertising or
    publicity relating to the software without the specific, prior written
    permission of Sam Leffler and Silicon Graphics.  

    THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND, 
    EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY 
    WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  

    IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR
    ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
    OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
    WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF 
    LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
    OF THIS SOFTWARE."

End-copyright-notice-for-Libtiff
*********************************************
*/

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <gtk/gtk.h>
#include <stdlib.h>
#include <string.h>
#include "interface.h"
#include "support.h"
#include "main.h"


#include <glib.h>
#include <glib/gprintf.h>
#include <glib/gstdio.h>
#include <assert.h>
#include <stdio.h>

#include <unistd.h>

#define gtkimageviewer 1
#define SECONDS_PER_DAY 86400
#define NUCLEUS_CHANNEL_TAG

#define MYDEBUG 1
#undef  MYDEBUG   

// Code from Andy Gordon CellID
#include "date_and_time.c"


enum{
  STRING_COL = 0,
  STRINGUTF8_COL,
  TIMETAG_COL,
  ISFILENAME_COL,
  NUM_COLS
} ;

enum{
 SORTID_NAME=0
};

PARAMS p;

GtkWidget *app1;
GtkWidget *inputDialog;
GtkWidget *segmentDialog;
GtkWidget *imageSetupDialog;
GtkWidget *setupDialog;
GtkWidget *calculatingDialog;
GtkWidget *aboutDialog;

GtkTreeStore      *treestore;
GtkTreeModel      *sortmodel;

GtkAnimView      *view  = NULL;
GtkIImageTool *dragger  = NULL;
GtkIImageTool *selector = NULL;
GtkIImageTool *painter  = NULL;


//item 0 is bright field image list
//subsequent items correspond to different fluorescence channel lists
GPtrArray *gpaImages = NULL; 

gboolean b_coming_from_multipos=FALSE;

int
main (int argc, char *argv[])
{
 GtkWidget *w;
 gint context_id;

 
#ifdef ENABLE_NLS
  bindtextdomain (GETTEXT_PACKAGE, PACKAGE_LOCALE_DIR);
  bind_textdomain_codeset (GETTEXT_PACKAGE, "UTF-8");
  textdomain (GETTEXT_PACKAGE);
#endif

  gtk_set_locale ();
  gtk_init (&argc, &argv);

  add_pixmap_directory (PACKAGE_DATA_DIR "/" PACKAGE "/pixmaps");

  app1             = create_app1 ();
  inputDialog      = create_inputDialog ();
  segmentDialog    = create_segmentDialog ();
  imageSetupDialog = create_ImageSetupDialog ();
  setupDialog      = create_setupDialog ();
  calculatingDialog= create_calculatingDialog();
  aboutDialog      = create_aboutdialog1();
    
  w          = lookup_widget(app1,"statusbar1");
  context_id = gtk_statusbar_get_context_id (GTK_STATUSBAR (w),"VCell");

  //Initialize parameter structure with default values.
  InitParams();
  
  //Disabling LoadParams in this version v0.1
  /*
  if(LoadParams(p.param_file->str)==0){
     StatusMessage(context_id, 
                    "No '.parameters.txt' file found for initialization");
  }
  */
  
  //Initialize forms accordingly
  Param2InputDialog();
  Param2ImageSetupDialog();
  Param2SegmentationDialog();
  Param2SetupDialog();

  //create treestore object
  treestore = gtk_tree_store_new(NUM_COLS,G_TYPE_STRING,G_TYPE_STRING,G_TYPE_DOUBLE,G_TYPE_BOOLEAN);
  #ifdef gtkimageviewer
   view = GTK_ANIM_VIEW (gtk_anim_view_new ());
   w = lookup_widget(app1,"vpaned2")    ;
   gtk_container_add(GTK_CONTAINER (w), GTK_WIDGET(view));
   gtk_widget_set_size_request(GTK_WIDGET(view),100,100);

  
   gtk_paned_set_position(GTK_PANED(w),300);
   w = lookup_widget(app1,"hpaned2");
   gtk_paned_set_position(GTK_PANED(w),200);
  #endif

  gtk_widget_show_all(app1);
  
  // hide stuff
  w = lookup_widget(app1,"vboxFilter");
  g_object_set(GTK_WIDGET(w),"visible",FALSE,NULL);
  w = lookup_widget(app1,"vboxView");
  g_object_set(GTK_WIDGET(w),"visible",FALSE,NULL);
 
  w = lookup_widget(imageSetupDialog,"label50");
  g_object_set(GTK_WIDGET(w),"visible",FALSE,NULL);
  w = lookup_widget(imageSetupDialog,"fcbOut");
  g_object_set(GTK_WIDGET(w),"visible",FALSE,NULL);

  gtk_main ();

  //Cleaning TMP stuff before leaving
  CleanTMPdir();
  
  return 0; 
}




void CleanTMPdir(void){
  gboolean valid;
  GString *saux = g_string_new(NULL);
  GError  *error = NULL;
  GDir     *gdir;
  gchar    *pfname, *pchar1=NULL;
  //remove TMP folder
  valid = g_file_test(p.tmp_path->str,G_FILE_TEST_EXISTS);
  valid = valid && g_file_test(p.tmp_path->str,G_FILE_TEST_IS_DIR);

  if(valid){
   gdir = g_dir_open(p.tmp_path->str,0,&error);
   while((pfname=g_dir_read_name(gdir))) {
    pchar1 = g_build_filename(p.tmp_path->str,pfname,NULL);
    valid = valid && (g_remove(pchar1)==0);
    g_free(pchar1);
   }
   g_dir_close(gdir); 
   valid = valid && (g_rmdir(p.tmp_path->str)==0);
   if(!valid){
       g_string_printf(saux,"Could not delete TMP folder:\n %s\nPlease remove it manually.",p.tmp_path->str);
       message_dialog(GTK_MESSAGE_WARNING,
                      GTK_WINDOW(app1),
                      "Warning",saux->str );
      g_string_free(saux,TRUE);
    }
  }  
}


void InitParams(void){
  gchar *pchar=NULL, 
        *pchar2=NULL;
  GRand *gr;
  guint32 i32;
  GString *saux = g_string_new(NULL);
  
  gr  = g_rand_new();
  i32 = g_random_int_range(1,10000);
  g_rand_free(gr);
  

  //create TMP dir to put files obtained during testing
  pchar = g_get_tmp_dir();
  g_string_printf(saux,"VCell%0d",i32);  
  pchar2= g_build_filename(pchar,saux->str,NULL);
  p.tmp_path = g_string_new(pchar2);
  g_free(pchar2);
  g_string_free(saux,TRUE);
  g_mkdir(p.tmp_path->str,0722);
  
  //of no use in this version v0.1
  p.param_file = g_string_new(NULL);
  g_string_printf(p.param_file,"%s",".parameters.txt");


  //R friendly output
  p.bPAW = FALSE;
 
  //of no use in this version v0.1 (user not allowed to change this)
  p.segment_values_t0               = FALSE; 
  p.max_dist_over_waist_t0          = -999.0;
  p.max_split_over_minor_t0         = -999.0;
  
  p.max_dist_over_waist             = 8.0;
  p.max_split_over_minor            = 0.5;
  
  p.min_pixels_per_cell               = 75;
  p.max_pixels_per_cell               = 1500;
  
  p.background_reject_factor          = 1.0;
  p.tracking_comparison               = 0.2;
  
  p.align_cells                       = CELL_ALIGN_CELL;
  
  p.image_type                        = BRIGHTFIELD;
  
  // v0.1 (user not allowed to change this)
  p.force_nucleus_in_center           = FALSE;
    
  p.splitted_fret_image               = FALSE;
  p.fret_nucleus                      = FRET_NUCLEUS_TOP;
  
  p.align_fl                          = ALIGN_FL_BF;
  
  // v0.1 (user not allowed to change this)
  p.bf_fl_mapping                     = BF_FL_MAPPING_FILENAME;

  p.treat_brightfield_as_fluorescence = 0;

  // v0.1 not implemented yet
  p.do_recomb  = 0;
  p.all_recomb = 0;
  p.nuc_recomb = 0;
  p.fl_recomb  = 0;

  if(p.bf_token==NULL)
    p.bf_token = g_string_new("BF");
  if(p.fl_token==NULL)
    p.fl_token = g_string_new("FP");  
  if(p.nucleus_channel==NULL)
    p.nucleus_channel = g_string_new(NULL);  
  if(p.flat_token==NULL)
    p.flat_token = g_string_new(NULL);  
  if(p.dark_token==NULL)
    p.dark_token = g_string_new(NULL);  
 
  p.bdark    = FALSE;
  p.bflat    = FALSE;
  p.bnucleus = FALSE; 
 
  p.same_path=TRUE;
  
  if(p.bf_path==NULL){
  #ifdef MYDEBUG
    pchar =g_get_current_dir();
    pchar2=g_build_filename(pchar,"images",NULL);
    p.bf_path = g_string_new(pchar2);
    g_free(pchar);
    g_free(pchar2);
   #else 
    p.bf_path = g_string_new(pchar=g_get_current_dir());
    g_free(pchar);
   #endif 
  }else{  
   g_string_assign(p.bf_path,pchar=g_get_current_dir());
   g_free(pchar);
  } 
  if(p.fl_path==NULL){
    p.fl_path = g_string_new(pchar=g_get_current_dir());
    g_free(pchar);
  }else{  
   g_string_assign(p.fl_path,pchar=g_get_current_dir());
   g_free(pchar);
  } 
  if(p.flat_path==NULL){
    p.flat_path = g_string_new(pchar=g_get_current_dir());
    g_free(pchar);
  }else{  
   g_string_assign(p.flat_path,pchar=g_get_current_dir());
   g_free(pchar);
  } 
  if(p.dark_path==NULL){
    p.dark_path = g_string_new(pchar=g_get_current_dir());
    g_free(pchar);
  }else{  
   g_string_assign(p.dark_path,pchar=g_get_current_dir());
   g_free(pchar);
  } 
  
  if(p.exec_file==NULL){
   #ifdef  G_OS_WIN32
     p.exec_file = g_string_new(strcat(g_path_get_dirname(argv0),"/cell.exe"));
   #else
    #ifdef MYDEBUG
     p.exec_file = g_string_new("/home/ariel/Projects/CellID/cell_id_1.4a/cell");
     if(!g_file_test(p.exec_file->str,G_FILE_TEST_EXISTS)){
       g_string_assign(p.exec_file,pchar=g_get_current_dir());
       g_free(pchar);
     }    
    #else
     p.exec_file = g_string_new(pchar=g_get_current_dir());
     g_free(pchar);
    #endif
   #endif
  }else{  
    g_string_assign(p.exec_file,pchar=g_get_current_dir());
    g_free(pchar);
  }     
  if(p.output_basename==NULL){
   p.output_basename = g_string_new("out");
  }else{  
   g_string_assign(p.output_basename,"out");
  }
  if(p.pos_token==NULL){
   p.pos_token = g_string_new("Position");
  }else{
   g_string_assign(p.pos_token,"Position");
  }
  if(p.time_token==NULL){
   p.time_token = g_string_new("time");
  }else{
   g_string_assign(p.time_token,"time");
  }
  if(p.separator==NULL){
   p.separator = g_string_new("_");
  }else{
   g_string_assign(p.separator,"_");
  }
  
  p.bpos =TRUE;
  p.btime=TRUE;
  p.bsep =TRUE;
  
  //will point to bf_path/pos_token
  p.output_path=g_string_new(NULL);
}

void SegmentationDialog2Param(void){
  GtkWidget    *w;
  
  w = lookup_widget(segmentDialog,"sbMaxSplit"); 
  p.max_split_over_minor = gtk_spin_button_get_value(
                              GTK_SPIN_BUTTON(w));
  w = lookup_widget(segmentDialog,"sbMaxDist"); 
  p.max_dist_over_waist = gtk_spin_button_get_value(
                              GTK_SPIN_BUTTON(w));
  
  w = lookup_widget(segmentDialog,"sbMinPix"); 
  p.min_pixels_per_cell = gtk_spin_button_get_value(
                              GTK_SPIN_BUTTON(w));
  w = lookup_widget(segmentDialog,"sbMaxPix"); 
  p.max_pixels_per_cell = gtk_spin_button_get_value(
                              GTK_SPIN_BUTTON(w));
  
  w = lookup_widget(segmentDialog,"sbBkgdRejectFactor"); 
  p.background_reject_factor = gtk_spin_button_get_value(
                              GTK_SPIN_BUTTON(w));
  w = lookup_widget(segmentDialog,"sbTracking"); 
  p.tracking_comparison = gtk_spin_button_get_value(
                              GTK_SPIN_BUTTON(w));
  
  p.align_cells = CELL_ALIGN_NO;
  w = lookup_widget(segmentDialog,"rbAlignCells"); 
  if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)) ) 
         p.align_cells = CELL_ALIGN_CELL;
  w = lookup_widget(segmentDialog,"rbAlignCellsBoundaries"); 
  if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)) ) 
         p.align_cells = CELL_ALIGN_BOUND;
         
  p.align_fl = ALIGN_NO;       
  w = lookup_widget(segmentDialog,"rbFirstAlignement"); 
  if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)) ) 
         p.align_fl = ALIGN_FL_t0;
  w = lookup_widget(segmentDialog,"rbBFAlignement");
  if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)) ) 
         p.align_fl = ALIGN_FL_BF;

}

void InputDialog2Param(void){
 GtkWidget    *w;
 gchar *saux;

 //saux=(gchar *)calloc(100,sizeof(gchar));
 ///////////////////
 // Image File tab//
 w = lookup_widget(inputDialog,"entryBF"); 
 g_string_assign(p.bf_token,gtk_entry_get_text(GTK_ENTRY(w)));
 w = lookup_widget(inputDialog,"entryFL"); 
 g_string_assign(p.fl_token,gtk_entry_get_text(GTK_ENTRY(w)));

 w = lookup_widget(inputDialog,"entryFlat"); 
 g_string_assign(p.flat_token,gtk_entry_get_text(GTK_ENTRY(w)));
 w = lookup_widget(inputDialog,"entryDark"); 
 g_string_assign(p.dark_token,gtk_entry_get_text(GTK_ENTRY(w)));
  
 w = lookup_widget(inputDialog,"fcbBrightField"); 
 saux = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(w));
 g_string_assign(p.bf_path,saux);
 g_free(saux);
 w = lookup_widget(inputDialog,"fcbFluo"); 
 saux=gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(w));
 g_string_assign(p.fl_path,saux);
 g_free(saux);
 w = lookup_widget(inputDialog,"fcbDark"); 
 saux=gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(w));
 g_string_assign(p.dark_path,saux);
 g_free(saux);     
 w = lookup_widget(inputDialog,"fcbFlat"); 
 saux=gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(w));
 g_string_assign(p.flat_path,saux);
 g_free(saux);               
 w = lookup_widget(inputDialog,"cbForcePath"); 
 p.same_path = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
 
 w = lookup_widget(inputDialog,"cbPosToken"); 
 p.bpos = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
 w = lookup_widget(inputDialog,"cbTimeToken"); 
 p.btime = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
 w = lookup_widget(inputDialog,"cbSeparator"); 
 p.bsep = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)); 
 w = lookup_widget(inputDialog,"entryPosToken"); 
 g_string_assign(p.pos_token,gtk_entry_get_text(GTK_ENTRY(w)));
 w = lookup_widget(inputDialog,"entryTimeToken"); 
 g_string_assign(p.time_token,gtk_entry_get_text(GTK_ENTRY(w)));
 w = lookup_widget(inputDialog,"entrySeparator"); 
 g_string_assign(p.separator,gtk_entry_get_text(GTK_ENTRY(w)));


 w = lookup_widget(inputDialog,"cboxMapping");
 p.bf_fl_mapping = GetComboBoxItem(GTK_COMBO_BOX(w),NULL);

 w          = lookup_widget(inputDialog,"cbBackground");
 p.bdark    = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
 w          = lookup_widget(inputDialog,"cbGradient");
 p.bflat    = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

 //free(saux);
 return;
}

void SetupDialog2Param(void){
 GtkWidget    *w;
 w = lookup_widget(setupDialog,"fcbExec");                 
 g_string_assign(p.exec_file,
                 gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(w)));
}

void ImageSetupDialog2Param(void){
 GtkWidget    *w;

 gint i;
 
 ////////////////////
 // image type tab //
 w = lookup_widget(imageSetupDialog,"cboxImageType");
 i = GetComboBoxItem(GTK_COMBO_BOX(w),NULL);
 p.image_type = i;
 
 w = lookup_widget(imageSetupDialog,"cbFret");
 if(!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))){
  p.splitted_fret_image=FALSE;
 }else{
  w = lookup_widget(imageSetupDialog,"cboxFretNucleus");
  p.fret_nucleus       = GetComboBoxItem(GTK_COMBO_BOX(w),NULL);
  p.splitted_fret_image=TRUE;
 }
 
 w = lookup_widget(imageSetupDialog,"cbBFasFL");
 p.treat_brightfield_as_fluorescence =
      gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

 w          = lookup_widget(imageSetupDialog,"cbNucleus");
 p.bnucleus = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

 w = lookup_widget(imageSetupDialog,"entryNucleus"); 
 g_string_assign(p.nucleus_channel,gtk_entry_get_text(GTK_ENTRY(w)));
 //////////////
 // Misc tab //
 w = lookup_widget(imageSetupDialog,"entryBasename");
 g_string_assign(p.output_basename,
                 gtk_entry_get_text(GTK_ENTRY(w)));

 w = lookup_widget(imageSetupDialog,"rbPAW");
 p.bPAW = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

 return;
}

void Param2InputDialog(){
 GtkWidget    *w;
 gboolean     b;

  
 /////////////
 // mapping // 
 w = lookup_widget(inputDialog,"cboxMapping");
 if(p.bf_fl_mapping==BF_FL_MAPPING_FILENAME){  
  SetComboBoxItem(GTK_COMBO_BOX(w),"file name pattern");
 }
 if(p.bf_fl_mapping==BF_FL_MAPPING_TIMETAG){  
  SetComboBoxItem(GTK_COMBO_BOX(w),"time (metamorph)");
 }
 if(p.bf_fl_mapping==BF_FL_MAPPING_LIST){  
  SetComboBoxItem(GTK_COMBO_BOX(w),"mapping file");
 }

 
 //////////////////////
 //   same path      //
 w = lookup_widget(inputDialog,"cbForcePath");
 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),
             p.same_path);
 
 w = lookup_widget(inputDialog,"fcbBrightField");
 gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER(w),
                                      (p.bf_path)->str );
 if(p.same_path){
   w = lookup_widget(inputDialog,"fcbFluo");
   gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER(w),
                                      (p.bf_path)->str );

   w = lookup_widget(inputDialog,"fcbFlat");
   gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER(w),
                                      (p.bf_path)->str ); 
   w = lookup_widget(inputDialog,"fcbDark");
   gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER(w),
                                      (p.bf_path)->str ); 
 }else{
   w = lookup_widget(inputDialog,"fcbFluo");
   gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER(w),
                                      (p.fl_path)->str );

   w = lookup_widget(inputDialog,"fcbFlat");
   gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER(w),
                                      (p.flat_path)->str ); 
   w = lookup_widget(inputDialog,"fcbDark");
   b = gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER(w),
                                      (p.dark_path)->str );  
 }
 
 /////////////
 // tokens  //
 w = lookup_widget(inputDialog,"entryBF");
 gtk_entry_set_text(GTK_ENTRY(w),(p.bf_token)->str);
 w = lookup_widget(inputDialog,"entryFL");
 gtk_entry_set_text(GTK_ENTRY(w),(p.fl_token)->str);

 gtk_widget_set_sensitive(GTK_WIDGET(w),p.bnucleus);
 w = lookup_widget(inputDialog,"entryFlat");
 gtk_entry_set_text(GTK_ENTRY(w),(p.flat_token)->str);
 gtk_widget_set_sensitive(GTK_WIDGET(w),p.bflat);
 w = lookup_widget(inputDialog,"entryDark");
 gtk_entry_set_text(GTK_ENTRY(w),(p.dark_token)->str);
 gtk_widget_set_sensitive(GTK_WIDGET(w),p.bdark);
 
 w = lookup_widget(inputDialog,"cbBackground");
 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),p.bdark);
 w = lookup_widget(inputDialog,"cbGradient");
 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),p.bflat);

 /////////////////
 // more tokens //
 w = lookup_widget(inputDialog,"entryPosToken");
 gtk_entry_set_text(GTK_ENTRY(w),(p.pos_token)->str);
 w = lookup_widget(inputDialog,"entryTimeToken");
 gtk_entry_set_text(GTK_ENTRY(w),(p.time_token)->str);
 w = lookup_widget(inputDialog,"entrySeparator");
 gtk_entry_set_text(GTK_ENTRY(w),(p.separator)->str);

 w = lookup_widget(inputDialog,"cbPosToken"); 
 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),p.bpos);
 w = lookup_widget(inputDialog,"cbTimeToken"); 
 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),p.btime);
 w = lookup_widget(inputDialog,"cbSeparator"); 
 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),p.bsep);  
 
}

void Param2SetupDialog(){
 GtkWidget    *w;
 w = lookup_widget(setupDialog,"fcbExec");
 if(g_file_test(p.exec_file->str,G_FILE_TEST_EXISTS))
  gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(w),
                                    (p.exec_file)->str);
}

void Param2ImageSetupDialog(){
 GtkWidget    *w;
 gchar        saux[100],saux2[100];
 gboolean     b;

 ////////////////
 // image type //
 w     = lookup_widget(imageSetupDialog,"cboxImageType");
 if(p.image_type==BRIGHTFIELD) sprintf(saux,"%s","brightfield");
 if(p.image_type==CONFOCAL_TRANSMISSION)sprintf(saux,"%s","confocal_transmission");
 if(p.image_type==METAMORPH_DECONVOLUTION)sprintf(saux,"%s","metamorph_deconvolution");
 if(SetComboBoxItem(GTK_COMBO_BOX(w), saux)==FALSE){
  //    message_dialog(GTK_MESSAGE_WARNING,
  //                   GTK_WINDOW(app1),
  //                   "Warning",
  //                   "No image type  for initialization\n");
 }
 
 //////////
 // FRET //
 w = lookup_widget(imageSetupDialog,"cbFret");
 if(p.splitted_fret_image==FALSE){  
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),FALSE);
 }else{
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),TRUE);
  if(p.fret_nucleus==FRET_NUCLEUS_TOP){
   sprintf(saux2,"%s","top");
  }
  if(p.fret_nucleus==FRET_NUCLEUS_BOTTOM){
   sprintf(saux2,"%s","bottom");
  }
  w = lookup_widget(imageSetupDialog,"cboxFretNucleus");
  SetComboBoxItem(GTK_COMBO_BOX(w), saux2);
 }  
  
 w = lookup_widget(imageSetupDialog,"cbNucleus");
 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),p.bnucleus);
 w = lookup_widget(imageSetupDialog,"entryNucleus");
 gtk_entry_set_text(GTK_ENTRY(w),(p.nucleus_channel)->str);
 
 //////////////
 // bf as fl //
 w = lookup_widget(imageSetupDialog,"cbBFasFL");
 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),
             p.treat_brightfield_as_fluorescence);
 
 w = lookup_widget(imageSetupDialog,"rbPAW");
 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),
             p.bPAW);



 //////////////
 // out path //
 w = lookup_widget(imageSetupDialog,"entryBasename");
 gtk_entry_set_text(GTK_ENTRY(w),(p.output_basename)->str); 
 
// w = lookup_widget(imageSetupDialog,"fcbOut");
// b=gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(w),
//                                        (p.output_path)->str);
 if(b==FALSE){
 
 }
 

}

void Param2SegmentationDialog(void){
 GtkWidget    *w;
 gchar        saux[100];

 ///////////////////////////////////
 // max split and dist over waist //
 w = lookup_widget(segmentDialog,"sbMaxSplit"); 
 gtk_spin_button_set_value(GTK_SPIN_BUTTON(w),p.max_split_over_minor);
 w = lookup_widget(segmentDialog,"sbMaxDist"); 
 gtk_spin_button_set_value(GTK_SPIN_BUTTON(w),p.max_dist_over_waist);
 
 //////////////////
 // pix per cell //
 w = lookup_widget(segmentDialog,"sbMinPix"); 
 gtk_spin_button_set_value(GTK_SPIN_BUTTON(w),p.min_pixels_per_cell);
 
 w = lookup_widget(segmentDialog,"sbMaxPix"); 
 gtk_spin_button_set_value(GTK_SPIN_BUTTON(w),p.max_pixels_per_cell);
 
 /////////////////
 // bkgd reject //
 w = lookup_widget(segmentDialog,"sbBkgdRejectFactor"); 
 gtk_spin_button_set_value(GTK_SPIN_BUTTON(w),p.background_reject_factor);
 
 //////////////
 // tracking //
 w = lookup_widget(segmentDialog,"sbTracking"); 
 gtk_spin_button_set_value(GTK_SPIN_BUTTON(w),p.tracking_comparison);
  
 //////////////////////
 //  cell alignment //
 if(p.align_cells==CELL_ALIGN_NO){
   w = lookup_widget(segmentDialog,"rbNoAlignCells"); 
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),TRUE);
 }  
 
 if(p.align_cells==CELL_ALIGN_CELL){
   w = lookup_widget(segmentDialog,"rbAlignCells"); 
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),TRUE);
 }  
 
 if(p.align_cells==CELL_ALIGN_BOUND){
   w = lookup_widget(segmentDialog,"rbAlignCellsBoundaries"); 
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),TRUE);
 }  
 
  
 ////////////////
 // Alignement //
 if(p.align_fl==ALIGN_NO){
   w = lookup_widget(segmentDialog,"rbNoAlignement"); 
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),TRUE);
 }  
 
 if(p.align_fl==ALIGN_FL_t0){
   w = lookup_widget(segmentDialog,"rbFirstAlignement"); 
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),TRUE);
 }  
 
 if(p.align_fl==ALIGN_FL_BF){
   w = lookup_widget(segmentDialog,"rbBFAlignement"); 
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),TRUE);
 }  

 ////////////
 // REcomb //
if(p.all_recomb){
   w = lookup_widget(segmentDialog,"cbNewRecomb"); 
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),TRUE);
 }
 if(p.nuc_recomb){
   w = lookup_widget(segmentDialog,"cbNucRecomb"); 
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),TRUE);
   
   w = lookup_widget(segmentDialog,"entryNucRecomb"); 
   sprintf(saux,"%f",p.nuc_recomb_val); 
   gtk_entry_set_text(GTK_ENTRY(w),saux);
   
   w = lookup_widget(segmentDialog,"cboxNucRecomb");  
 //FALTA YA DEBERIAN ESTAR DETERMINADOS LOS TIPOS
 //DE FLUORESCENCIA y CON esos se deberia llenar esto
 //como saberlo desde el parameter.txt?
 }

 if(p.fl_recomb){
   w = lookup_widget(segmentDialog,"cbFLRecomb"); 
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),TRUE);
   
   w = lookup_widget(segmentDialog,"entryFLRecomb"); 
   sprintf(saux,"%f",p.fl_recomb_val); 
   gtk_entry_set_text(GTK_ENTRY(w),saux);
   
   w = lookup_widget(segmentDialog,"cboxFLRecomb");  
   //FALTA
 }

}

gboolean LoadParams(gchar *fname){

   FILE *fp;
   char line[500],line2[500],line3[500];
   double double_tmp;
   int    itmp;

   fp=fopen(p.param_file->str,"r");
   if(fp==NULL) return(0);

   p.align_cells      = CELL_ALIGN_NO;
   p.align_fl         = ALIGN_NO;
   p.splitted_fret_image = FALSE;
   g_string_assign(p.nucleus_channel,"");
   p.do_recomb        = FALSE;
   p.nuc_recomb       = FALSE;
   p.all_recomb       = FALSE;
   p.fl_recomb        = FALSE;
   p.force_nucleus_in_center = FALSE;
   p.segment_values_t0   = FALSE;
   while((fgets(line,450,fp))!=NULL){ //next line (while not EOF)
      if (line[0]!='#'){ //if not a comment
	if ((strstr(line,"max_split_over_minor"))!=NULL){
          if ((strstr(line,"_t0"))!=NULL){
	          sscanf(line,"%s %le",line2,&double_tmp);
	          p.max_split_over_minor_t0=double_tmp;
	          p.segment_values_t0=TRUE;
       	  }else{
	          sscanf(line,"%s %le",line2,&double_tmp);
	          p.max_split_over_minor=double_tmp;
	  }
	}else if((strstr(line,"max_dist_over_waist"))!=NULL){
	        if ((strstr(line,"_t0"))!=NULL){
	          sscanf(line,"%s %le",line2,&double_tmp);
	          p.max_dist_over_waist_t0=double_tmp;
	        }else{
	           p.segment_values_t0=TRUE;
	           sscanf(line,"%s %le",line2,&double_tmp);
	           p.max_dist_over_waist=double_tmp;
	        }
	}else if((strstr(line,"background_reject_factor"))!=NULL){
	         sscanf(line,"%s %le",line2,&double_tmp);
	         p.background_reject_factor=double_tmp;
	}else if((strstr(line,"tracking_comparison"))!=NULL){
	         sscanf(line,"%s %le",line2,&double_tmp);
	         p.tracking_comparison=double_tmp;
	}else if((strstr(line,"max_pixels_per_cell"))!=NULL){
	         sscanf(line,"%s %i",line2,&itmp);
	         p.max_pixels_per_cell=itmp;
	}else if((strstr(line,"min_pixels_per_cell"))!=NULL){
	         sscanf(line,"%s %i",line2,&itmp);
	         p.min_pixels_per_cell=itmp;
	}else if((strstr(line,"force_nucleus_in_center"))!=NULL){
	         p.force_nucleus_in_center=TRUE;
	}else if(strstr(line,"image_type")!=NULL){
	         if(strstr(line,"brightfield")!=NULL)
	           p.image_type = BRIGHTFIELD;
	         if(strstr(line,"metamorph")!=NULL) 
	           p.image_type = METAMORPH_DECONVOLUTION;
	         if(strstr(line,"confocal")!=NULL)
	           p.image_type = CONFOCAL_TRANSMISSION;
        //V1.2a bright field to fluorescence mapping
        }else if(strstr(line,"bf_fl_mapping")!=NULL){
                sscanf(line,"%s %s",line2,line3);
                if(strstr(line3,"time")!=NULL){
                  //printf("Mapping bright field and fluorescence images by time.\n");
                  p.bf_fl_mapping=BF_FL_MAPPING_TIMETAG;
                } else if((strstr(line3,"list")!=NULL)||(strstr(line3,"file_name")!=NULL)){
                  //printf("Mapping bright field and fluorescence images by list order.\n");
                  //printf("Same number of elemtes requiered in list files.\n");
                  p.bf_fl_mapping=BF_FL_MAPPING_FILENAME;
                } else {
                  /*TODO 
                  IssueWarn
                  printf("-%s- is a invalid value for bf_fl_mapping in parameter.txt.\n",line3);
                  printf("Using time mapping by default.\n");*/
                  p.bf_fl_mapping=BF_FL_MAPPING_TIMETAG; 
                }
	}else if(strstr(line,"fret")!=NULL){
	        sscanf(line,"%s %s",line2,line3);

	        //Check if should use top or bottom of image that's used
	        //for nuclear label. (If this isn't set, then it defaults
	        //to bottom of image if there is a third image set and to
	        //the top of the image otherwise--this is done below all this).
	        if(strstr(line3,"nuclear_top")!=NULL){
	          p.fret_nucleus=FRET_NUCLEUS_TOP;
	        }else if(strstr(line3,"nuclear_bottom")!=NULL){
	          p.fret_nucleus=FRET_NUCLEUS_BOTTOM;
	        }
                p.splitted_fret_image = TRUE;
	}else if(strstr(line,"align_fl_to_first")!=NULL){
	        //Check whether we want to use the first fluorescence image to
	        //align all the others
	        //printf("Will align all FL files to first FL file.\n");
	        p.align_fl=ALIGN_FL_t0;
	}else if(strstr(line,"align_fl_to_bf")!=NULL){
	        //or whether to align fl images to bf
	        //printf("Will align first FL files to brightfield.\n");
	        p.align_fl=ALIGN_FL_BF;
	        //Check for third list of images. The use of the
	        //images depends on the image_type.
	}else if(strstr(line,"treat_brightfield_as_fluorescence_also")!=NULL){
	        p.treat_brightfield_as_fluorescence=1;
	        //printf("Adding BF image as additional fluorescence image.\n");
        }else if (strstr(line,"nucleus_channel")!=NULL){
                sscanf(line,"%s %s",line2,line3);
	          g_string_assign(p.nucleus_channel,line3);
        }else if(strstr(line,"align_individual_cells")!=NULL){
	        p.align_cells=CELL_ALIGN_CELL;
	        //see if "boundary" is part of name
	        if(strstr(line,"align_individual_cells_boundary")!=NULL){
	          p.align_cells=CELL_ALIGN_BOUND;
	          //printf("Will wiggle each cell around to re-align with BF");
	          //printf(" using boundary.\n");
	        }//else{
	          //printf("Will wiggle each cell around to re-align with BF.\n");
	        //}
         }
       }  //End of check that first character wasn't a "#"
     } //End of while loop over parameters.txt
     fclose(fp);
     return(1);
}

gint StatusMessage(gint cid,gchar *msg){
 GtkWidget *w;
 
 w     = lookup_widget(app1,"statusbar1");
 gtk_statusbar_pop (GTK_STATUSBAR (w), cid);
 
 return(gtk_statusbar_push (GTK_STATUSBAR (w), cid, msg));
}

GtkResponseType message_dialog(GtkMessageType type,GtkWindow *parent,const gchar *title, const gchar *message){
        /*
        This function is used for all 4 types of GtkMessageDialog. See the API
        documentation for GtkMessageDialog for further information.

                type - 1 of 4 enum values for the type of message dialog.
                parent - The parent window for the dialog (window1)
                title - A string for the dialog box title
                message - A string to display in the dialog box
        */

        GtkWidget *dialog;
        GtkResponseType response;
        GtkButtonsType buttons;

        if (type == GTK_MESSAGE_QUESTION) buttons = GTK_BUTTONS_YES_NO;
        else buttons = GTK_BUTTONS_OK;

        /* create the dialog */
        dialog = gtk_message_dialog_new (parent,
                        GTK_DIALOG_DESTROY_WITH_PARENT,
                        type,
                        buttons,
                        message);

        /* set title */
        gtk_window_set_title(GTK_WINDOW(dialog), title);

        /* show dialog */
        response = gtk_dialog_run (GTK_DIALOG (dialog));
        gtk_widget_destroy (dialog);

        return response;
}

gboolean SetComboBoxItem(GtkComboBox *w,gchar *saux){
 GtkTreeModel *gtm;
 GtkTreeIter  iter;
 gint nitems,i; 
 gchar *sitem;
 gboolean a;
 
 gtm   = gtk_combo_box_get_model(GTK_COMBO_BOX(w));
 g_assert (gtm);

 nitems = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(gtm),NULL);
 a=gtk_tree_model_get_iter_first(GTK_TREE_MODEL(gtm),&iter);
 i=0;
 while(i<nitems &&            
       gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(gtm),&iter,NULL,i)){
  gtk_tree_model_get (gtm, &iter, 0, &sitem, -1);
  if(strcmp(sitem,saux)==0){
   break;
  }else i++;
 }
 if(i<nitems){
  gtk_combo_box_set_active(w, i);
  return(TRUE);
 }else{
  return(FALSE);
 }
}

gint GetComboBoxItem(GtkComboBox *w,gchar **saux){
 if(saux!=NULL)
   *saux = gtk_combo_box_get_active_text(GTK_COMBO_BOX(w));
 return(gtk_combo_box_get_active(GTK_COMBO_BOX(w)));
 /*
 GtkTreeModel *gtm;
 GtkTreeIter  iter;
 gint i=0,iactive;
 
 gtm     = gtk_combo_box_get_model(GTK_COMBO_BOX(w));
 iactive = gtk_combo_box_get_active(GTK_COMBO_BOX(w));
 gtk_tree_model_get_iter_first(gtm, &iter);
 while(iactive<i) gtk_tree_model_iter_next(gtm,&iter);
 */
}

void UpdateEntryFileName(void){
 GString *saux,*ssep;
 GtkWidget *w;
 
 saux = g_string_new("?");
 ssep = g_string_new("");
 w = lookup_widget(inputDialog,"cbSeparator");
 if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))){
  w = lookup_widget(inputDialog,"entrySeparator");
  g_string_assign(ssep,gtk_entry_get_text(GTK_ENTRY(w)));
 } 
 
 w = lookup_widget(inputDialog,"entryFL");
 g_string_append(saux,gtk_entry_get_text(GTK_ENTRY(w)));
 
 w = lookup_widget(inputDialog,"cbPosToken");
 if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))){
  w = lookup_widget(inputDialog,"entryPosToken");
  g_string_append(saux,ssep->str);
  g_string_append(saux,gtk_entry_get_text(GTK_ENTRY(w)));
  g_string_append(saux,"*");
 } 

 w = lookup_widget(inputDialog,"cbTimeToken");
 if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))){
  w = lookup_widget(inputDialog,"entryTimeToken");
  g_string_append(saux,ssep->str);
  g_string_append(saux,gtk_entry_get_text(GTK_ENTRY(w)));
  g_string_append(saux,"*");
 }
 g_string_append(saux,".tif");
   
 w = lookup_widget(inputDialog,"entryFileName");

 gtk_entry_set_text(GTK_ENTRY(w),saux->str);
 return;
}

GString *GetFileNamePattern(gchar *imgType){
 GString *saux,
         *ssep;
 GtkWidget *w;
 
 saux = g_string_new(imgType);
 ssep = g_string_new(NULL);
 w = lookup_widget(inputDialog,"cbSeparator");
 if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))){
  w = lookup_widget(inputDialog,"entrySeparator");
  g_string_assign(ssep,gtk_entry_get_text(GTK_ENTRY(w)));
 } 
 
 w = lookup_widget(inputDialog,"cbPosToken");
 if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))){
  w = lookup_widget(inputDialog,"entryPosToken");
  g_string_append(saux,ssep->str);
  g_string_append(saux,gtk_entry_get_text(GTK_ENTRY(w)));
  g_string_append(saux,"*");
 } 

 w = lookup_widget(inputDialog,"cbTimeToken");
 if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))){
  w = lookup_widget(inputDialog,"entryTimeToken");
  g_string_append(saux,ssep->str);
  g_string_append(saux,gtk_entry_get_text(GTK_ENTRY(w)));
  g_string_append(saux,"*");
 }
 g_string_append(saux,".tif");
   
 g_string_free(ssep,TRUE);  
 return(saux);
}

gboolean CheckSegmentationForm(void){
 return(TRUE);
}

gboolean CheckInputForm(void){
 GtkWidget *w;
 gchar *pchar;
 gchar *saux;
 gboolean valid;
 
 valid = TRUE;
 w = lookup_widget(inputDialog,"cbGradient");
 if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))){
   w = lookup_widget(inputDialog,"entryFlat");
   pchar = gtk_entry_get_text(GTK_ENTRY(w));
   if(pchar[0]=='\0'){
    valid = FALSE;
    message_dialog(GTK_MESSAGE_WARNING,
                     GTK_WINDOW(app1),
                     "Warning",
                     "Uncheck flat file option or set proper basename");
   }
 }
 w = lookup_widget(inputDialog,"cbBackground");
 if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))){
   w = lookup_widget(inputDialog,"entryDark");
   pchar = gtk_entry_get_text(GTK_ENTRY(w));
   if(pchar[0]=='\0'){
    valid = FALSE;
    message_dialog(GTK_MESSAGE_WARNING,
                     GTK_WINDOW(app1),
                     "Warning",
                     "Uncheck dark file option or set proper basename");
   }
 }

 return(valid);
}

gboolean CheckImageSetupForm(void){
  gboolean valid=TRUE;
  GtkWidget *w;
  GString *saux=g_string_new(NULL);

  w = lookup_widget(imageSetupDialog,"entryNucleus");
  g_string_assign(saux,gtk_entry_get_text(GTK_ENTRY(w)));
 
 //check 3 letter code. First one channel - last two fl_token
 valid = valid && 
         (saux->len==3) && (strcmp(&(saux->str)[1],p.fl_token->str)==0);
 
  return(valid);
};

gboolean CheckExecFile(void){
  FILE *fp=NULL;
  gchar *pchar=NULL;
  GString *saux=g_string_new(NULL);
  GtkWidget *w;
  
  w = lookup_widget(setupDialog,"fcbExec");
  pchar=gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(w));
  
  if(!g_file_test(pchar,G_FILE_TEST_IS_REGULAR)){
      g_string_printf(saux,"%s executable file could not be found",pchar);
      message_dialog(GTK_MESSAGE_ERROR,
                     GTK_WINDOW(app1),
                     "Error",
                     saux->str);
      return(FALSE);               
  }
  return(TRUE);
  /*
  if(pchar) fp=fopen(pchar,"r");
  if((fp==NULL) || (pchar==NULL)){
      g_string_printf(saux,"%s executable file could not be found",pchar);
      message_dialog(GTK_MESSAGE_ERROR,
                     GTK_WINDOW(app1),
                     "Error",
                     saux->str);
      return(FALSE);               
  }
  if(fp) fclose(fp);
  return(TRUE);
  */
};


gint SortFileNames__A(gconstpointer a, gconstpointer b){

 IMG *pa,*pb;
 gint i;
 
 pa=(IMG *)a;
 pb=(IMG *)b;
 
 //pos sorting
 i = strcmp(pa->pos,pb->pos);
 if(i==0){
  if(pa->time){
   return( (strcmp(pa->time,pb->time)>0)?1:-1);
  }else{ 
   return(0);
  } 
 }else{
  return(i);
 }
}




gint SortFileNames(gconstpointer a, gconstpointer b){

 IMG *pa,*pb;

 
 pa=(IMG *)a;
 pb=(IMG *)b;

 return(strcmp(pa->fname_utf,pb->fname_utf));
}


gint SortFileNamesTimeTag(gconstpointer a, gconstpointer b){

 IMG *pa,*pb;
 gint i;
 
 pa=(IMG *)a;
 pb=(IMG *)b;
 
 //pos sorting
 i = strcmp(pa->pos,pb->pos);
 if(i==0){
  if(pa->timetag){
   return( (pa->timetag > pb->timetag)?1:-1);
  }else{ 
   return(0);
  } 
 }else{
  return(i);
 }
}


//Funcion para renderear la info que 
// se muestra
void CellDataFunc(GtkTreeViewColumn *col,
                    GtkCellRenderer   *renderer,
                    GtkTreeModel      *model,
                    GtkTreeIter       *iter,
                    gpointer           user_data){
                    
  gchar *pchar;
  
  gtk_tree_model_get(model,iter,STRING_COL, &pchar,-1);
  g_object_set(renderer,"text",pchar,NULL);
  g_free(pchar); 
}

void FreeFileNameList(gpointer data,gpointer user_data){
 IMG *p = (IMG *)data;
 g_free(p->fname);
 g_free(p->pos);
 g_free(p->time);
 g_free(p->fname_utf);
 g_free(p->pos_utf);
 g_free(p->time_utf);

}

int GetTimeTag(GList *limg,gchar *spath){
 GList *n;
 int dtmp,it;
 float xstage, ystage;
 IMG *pimg;
 GString *saux=g_string_new(NULL);
 int bOK=1;
 gchar *s;
 double t;
 
 n = g_list_first(limg);
 while(n){
  pimg = (IMG *)n->data;
  s = g_build_filename(spath,pimg->fname,NULL);
  if(get_date_and_time(s,&dtmp,&it,&xstage,&ystage)==0){
     g_string_printf(saux,"Couldn't get date and time for %s.\n"
                          "Try switching to BF_FL_MAPPING_FILENAME mode\n",s);
     message_dialog(GTK_MESSAGE_ERROR,GTK_WINDOW(app1),
                  "Error",saux->str); 
     dtmp=0;
     t=0.0;
     xstage=-99999.0;
     ystage=-99999.0;
     bOK=0;
     g_free(s);     
     break;
  }
  t=((double)(it)/1000)/((float)SECONDS_PER_DAY) + (double)dtmp; 
  pimg->timetag = t ;
  n = n->next;
  g_free(s);
 }
 
 g_string_free(saux,TRUE);
 return(bOK);
}

// 1. Given images path and tokens, makes separate lists of 
//    IMG elements for brightfield, and fluorescence image filenames
// 2. Sort lists according timeTag
// 3. DARK y FLAT????
// 4. Build or Rebuild treestore
// 5. Connect treview <->treestore
// Note: sortmodel not implemented
gboolean LoadImageNames(void){
 GtkWidget  *w;
 GList      *bfFileNames   = NULL, 
            *flFileNames   = NULL,
            *flatFileNames = NULL,
            *darkFileNames = NULL,
            *thirdFileNames= NULL;
 GList      *n;

 GString *saux =NULL, 
         *saux2=g_string_new(NULL);
 gint iAlready;
 IMG *pimg;     
  
 GtkTreeViewColumn   *col;
 GtkCellRenderer     *renderer;
 void (*CreateAndFillModel)(GList *,gboolean);
 gint (*SortFileNamesFunc)(gconstpointer a, gconstpointer b);
 
 //bright field
 saux        = GetFileNamePattern(p.bf_token->str);
 bfFileNames = GetFileNames(p.bf_path->str,saux->str);
 g_string_free(saux,TRUE);
 SortFileNamesFunc = SortFileNames; //sort by name by default
 if(p.bf_fl_mapping==BF_FL_MAPPING_TIMETAG){
   if(GetTimeTag(bfFileNames,p.bf_path->str)==0){
     return(FALSE);
   }else{
    SortFileNamesFunc = SortFileNamesTimeTag;
   }
 } 
 bfFileNames=g_list_sort(bfFileNames,SortFileNamesFunc);
 
 #ifdef MYDEBUG  
 gchar *pchar1;
 n = g_list_first(bfFileNames);
 while(n){
  pchar1 = ((IMG *)(n->data))->fname;
  g_printf("%d-%s %f\n",g_list_length(bfFileNames),pchar1,((IMG *)(n->data))->timetag);
  n=n->next;
 }
 #endif
 
 //fluor
 g_string_printf(saux2,"?%s",p.fl_token->str);
 saux  = GetFileNamePattern(saux2->str);
 w = lookup_widget(inputDialog,"cbForcePath");
 if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)))
  g_string_assign(saux2,p.bf_path->str);
 else
  g_string_assign(saux2,p.fl_path->str); 
 flFileNames = GetFileNames(saux2->str,saux->str);
 if(p.bf_fl_mapping==BF_FL_MAPPING_TIMETAG){
   if(GetTimeTag(flFileNames,saux2->str)==0){
     return(FALSE);
   }else{
    SortFileNamesFunc = SortFileNamesTimeTag;
   }
 } 
 flFileNames=g_list_sort(flFileNames,SortFileNamesFunc);
 
 #ifdef MYDEBUG  
 n = g_list_first(flFileNames);
 while(n){
  pchar1 = ((IMG *)(n->data))->fname;
  g_printf("%s %f\n",pchar1,((IMG *)(n->data))->timetag);
  n=n->next;
 }
 #endif
 

 // dark y flat files
 if(p.bdark){
  g_string_printf(saux,"*%s*",p.dark_token->str);
  w = lookup_widget(inputDialog,"cbForcePath");
  if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)))
   g_string_assign(saux2,p.bf_path->str);
  else
   g_string_assign(saux2,p.dark_path->str); 
  darkFileNames = GetFileNames(saux2->str,saux->str);
  n = g_list_first(darkFileNames);
  if(n==NULL){
    g_string_prepend(saux,"No dark file image match ");
    g_string_append(saux ,"\nDisabling dark file option.");
    message_dialog(GTK_MESSAGE_WARNING,GTK_WINDOW(app1),
                  "Warning",saux->str);  
    p.bdark = FALSE;
    w = lookup_widget(inputDialog,"cbBackground");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),FALSE);
  }
  #ifdef MYDEBUG  
  while(n){
   pchar1 = ((IMG *)(n->data))->fname;
   g_printf("%s\n",pchar1);
   n=n->next;
  }  
  #endif
 }
 
 if(p.bflat){
  g_string_printf(saux,"*%s*",p.flat_token->str);
  w = lookup_widget(inputDialog,"cbForcePath");
  if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)))
   g_string_assign(saux2,p.bf_path->str);
  else
   g_string_assign(saux2,p.flat_path->str); 
  flatFileNames = GetFileNames(saux2->str,saux->str);
  n = g_list_first(flatFileNames);
  if(n==NULL){
    g_string_prepend(saux,"No flat file images match ");
    g_string_append(saux ,"\nDisabling flat file option.");
    message_dialog(GTK_MESSAGE_WARNING,GTK_WINDOW(app1),
                  "Warning",saux->str);  
    p.bflat = FALSE;
    w = lookup_widget(inputDialog,"cbGradient");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w),FALSE);
  }  
  #ifdef MYDEBUG  
  while(n){
   pchar1 = ((IMG *)(n->data))->fname;
   g_printf("%s\n",pchar1);
   n=n->next;
  }  
  #endif
 }
 
 g_string_free(saux,TRUE);
 g_string_free(saux2,TRUE);
 
 //already been here?
 iAlready = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(treestore),NULL);
 if(iAlready){
   //clear datastore contents
   gtk_tree_store_clear(treestore);
   
   //make sensitive image-related tool buttons
   gtk_anim_view_set_anim (view, NULL);
   w = lookup_widget(app1,"tbZoomIn");
   gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);
   w = lookup_widget(app1,"tbZoomOut");
   gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);
   w = lookup_widget(app1,"tbZoom1");
   gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);
   w = lookup_widget(app1,"tbZoomFit");
   gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);
   
 }
 
 
 //Create and fill data store: pos + time || pos || time || 0
 n=g_list_first(bfFileNames);
 if(n==NULL){
   message_dialog(GTK_MESSAGE_WARNING,GTK_WINDOW(app1),
                  "Error","No BF data to create data store");
   return(FALSE);               
 }else{
  pimg=((IMG *)(n->data));
  if(pimg->pos!=NULL && pimg->time!=NULL)
   CreateAndFillModel = CreateAndFillModel2;
                       
  if((pimg->pos!=NULL && pimg->time==NULL) ||
     (pimg->pos==NULL && pimg->time!=NULL) )
   CreateAndFillModel = CreateAndFillModel1;
                       /*
  if(pimg->pos==NULL && pimg->time==NULL)
   CreateAndFillModel0(bfFileNames, flFileNames,
                       thirdFileNames, darkFileNames,
                       flatFileNames); */
  } 
  CreateAndFillModel(bfFileNames,1);
  CreateAndFillModel(flFileNames,0);
  CreateAndFillModelX(darkFileNames,"Backgd Illumination"); //p.dark_token->str);
  CreateAndFillModelX(flatFileNames,"Uneven Illumination");//p.flat_token->str);
  

 if(!iAlready){
  sortmodel = gtk_tree_model_sort_new_with_model(GTK_TREE_MODEL(treestore));

  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(sortmodel), 
                                  STRINGUTF8_COL,
                                  sort_func_filename, 
                                  NULL, 
                                  NULL);

  gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(sortmodel),
                                         STRINGUTF8_COL, 
                                         GTK_SORT_ASCENDING);

 
  //Connect TreeDataStore with TreeView
  w = lookup_widget(app1,"treeview1");
  if(TRUE){
   gtk_tree_view_set_model(GTK_TREE_VIEW(w), GTK_TREE_MODEL(treestore));
   g_object_unref(sortmodel);
  }else{ 
   gtk_tree_view_set_model(GTK_TREE_VIEW(w), GTK_TREE_MODEL(sortmodel));
   g_object_unref(sortmodel);
  }
   
  g_signal_connect(w, "button-press-event", (GCallback) view_onButtonPressed, NULL);
  g_signal_connect(w, "popup-menu", (GCallback) view_onPopupMenu, NULL);
  
  
  //Create column 1
  col = gtk_tree_view_column_new();
  //gtk_tree_view_column_set_title(col, "First Name");
  /* pack tree view column into tree view */
  gtk_tree_view_append_column(GTK_TREE_VIEW(w), col);
  renderer = gtk_cell_renderer_text_new();
  /* pack cell renderer into tree view column */
  gtk_tree_view_column_pack_start(col, renderer, TRUE);
  /* connect a cell data function */
  gtk_tree_view_column_set_cell_data_func(col, renderer, CellDataFunc, NULL, NULL);
 }
 
 g_list_foreach(bfFileNames,FreeFileNameList,NULL);
 g_list_free(bfFileNames);
 g_list_foreach(flFileNames,FreeFileNameList,NULL);
 g_list_free(flFileNames);
 g_list_foreach(thirdFileNames,FreeFileNameList,NULL);
 g_list_free(thirdFileNames);
 g_list_foreach(flatFileNames,FreeFileNameList,NULL);
 g_list_free(flatFileNames);
 g_list_foreach(darkFileNames,FreeFileNameList,NULL);
 g_list_free(darkFileNames);
  
 return(TRUE);
}

void CreateAndFillModelX(GList *bf, gchar *ftoken){
 GtkTreeIter  iter1, iter2;
 GList *n;
 IMG   *pimg;
 gchar *pchar;
 
 n = bf;                                                         
 if(n==NULL) return;
 
 gtk_tree_model_get_iter_first(GTK_TREE_MODEL(treestore), &iter1);
 gtk_tree_store_insert_before(treestore,&iter1,NULL,NULL); 
 gtk_tree_store_set(treestore,&iter1,STRING_COL,ftoken,-1);
 pchar = g_utf8_collate_key (ftoken,-1);
 gtk_tree_store_set(treestore,&iter1,STRINGUTF8_COL, pchar,-1);
 g_free(pchar);
 gtk_tree_store_set(treestore,&iter1,ISFILENAME_COL, FALSE,-1);
 

 while(n){
  pimg=((IMG *)(n->data));
  gtk_tree_store_append(treestore, &iter2, &iter1);
  gtk_tree_store_set(treestore, &iter2,STRING_COL, pimg->fname,-1);
  pchar = g_utf8_collate_key (pimg->fname,-1);
  gtk_tree_store_set(treestore,&iter2,STRINGUTF8_COL, pchar,-1);
  g_free(pchar);
  gtk_tree_store_set(treestore,&iter2,ISFILENAME_COL, FALSE,-1);
  n = n->next;
 }
}

void CreateAndFillModel1(GList *bf, gboolean bCreate){
 GtkTreeIter  iter0, iter1, iter2;
 GList *n;
 IMG   *pimg;
 gchar *pchar,*pchar2;
 
 if(bCreate){
  gtk_tree_store_append(treestore, &iter0, NULL);
  gtk_tree_store_set(treestore,&iter0,STRING_COL , "root" ,-1);
 }else{
  gtk_tree_model_get_iter_first(GTK_TREE_MODEL(treestore),&iter0);
 }
 
 n = bf;
 while(n){
  pimg=((IMG *)(n->data));
  if(!GetTreeIter(pimg->pos, &iter1, &iter0)){ //if no level1 (pos), create it.
   gtk_tree_store_append(treestore, &iter1, &iter0);
   gtk_tree_store_set(treestore, &iter1,STRING_COL, pimg->pos,-1);
   pchar = g_utf8_collate_key (pimg->pos,-1);
   gtk_tree_store_set(treestore,&iter1,STRINGUTF8_COL, pchar,-1);
   g_free(pchar);
  }

  gtk_tree_store_append(treestore, &iter2, &iter1);
  gtk_tree_store_set(treestore,&iter2,STRING_COL , pimg->fname,-1);
  gtk_tree_store_set(treestore,&iter2,TIMETAG_COL, pimg->timetag,-1);
  pchar2= g_utf8_casefold(pimg->fname,-1);
  pchar = g_utf8_collate_key_for_filename (pchar2,-1);
  gtk_tree_store_set(treestore,&iter2,STRINGUTF8_COL, pchar,-1);
  g_free(pchar);
  g_free(pchar2);
  n = n->next;
 }
}

double GetIMGTimeTag(GList *n,gchar *fname){
 IMG *pimg;
  
 while(n){
  pimg=((IMG *)(n->data));
  if(strcmp(fname,pimg->fname)==0){
   return(pimg->timetag);
  }
  n = n->next;
 }
 return(-1);
}

// Build the 2 level treestore (root->posX->timeX)
// bCreate should be TRUE only for BF images.
void CreateAndFillModel2(GList *bf, gboolean bCreate){
 GtkTreeIter  iter0, iter1, iter2, iter3;
 GList *n=NULL;
 IMG   *pimg=NULL;
 gchar *pchar=NULL, *pchar2=NULL;
 int iGet;
 GString *saux = g_string_new(NULL);  
 int bout;
 double ttag;
 
 if(bCreate){
  gtk_tree_store_append(treestore, &iter0, NULL);
  gtk_tree_store_set(treestore,&iter0,STRING_COL , "root" ,-1);
 }else{
  gtk_tree_model_get_iter_first(GTK_TREE_MODEL(treestore),&iter0);
 }
 
 
 n = bf;
 while(n){
  pimg=((IMG *)(n->data));
  if(!GetTreeIter(pimg->pos, &iter1, &iter0)){ //if no level1 (pos), create it.
   gtk_tree_store_append(treestore, &iter1, &iter0);
   gtk_tree_store_set(treestore, &iter1,STRING_COL, pimg->pos,-1);
   pchar = g_utf8_collate_key (pimg->pos,-1);
   gtk_tree_store_set(treestore,&iter1,STRINGUTF8_COL, pchar,-1);
   g_free(pchar);
   gtk_tree_store_set(treestore,&iter1,ISFILENAME_COL, FALSE,-1);
  }
  if(p.bf_fl_mapping==BF_FL_MAPPING_FILENAME){
   iGet = GetClosestTimeTreeIter(pimg->time, &iter1, &iter2);
   if(iGet==0){
    if(bCreate){
     gtk_tree_store_append(treestore, &iter2, &iter1);
     //g_string_printf(saux,"%s %s",p.bf_token->str,pimg->time);
     gtk_tree_store_set(treestore   , &iter2,STRING_COL, pimg->time,-1);
     pchar = g_utf8_collate_key (pimg->time,-1);
     gtk_tree_store_set(treestore,&iter2,STRINGUTF8_COL, pchar,-1);
     g_free(pchar);
     gtk_tree_store_set(treestore,&iter2,ISFILENAME_COL, FALSE,-1);
    } 
   }else{
    #ifdef MYDEBUG
    if(iGet<0){
     //the BF list is supossed to be sorted so you shouldn't
     //be here for BF images
     g_string_printf(saux,"No previous BF image for %s",pimg->fname);
     message_dialog(GTK_MESSAGE_WARNING,GTK_WINDOW(app1),
                  "Warning",saux->str);
    }
    #endif 
   } 
  }else{ 
   //p.bf_fl_mapping==BF_FL_MAPPING_TIMETAG
   //different dataflow than for MAPPING_FILENAME
   //(better?)
    bout = (g_str_has_suffix(pimg->fname,".out.tif") ||
            g_str_has_suffix(pimg->fname,".out.TIF"));
    if(bCreate && !bout){ //only BF images
      gtk_tree_store_append(treestore, &iter2, &iter1);
      
      //g_string_printf(saux,"%s %s",p.bf_token->str,pimg->time);
      gtk_tree_store_set(treestore   , &iter2,STRING_COL, pimg->time,-1);
      gtk_tree_store_set(treestore   , &iter2,TIMETAG_COL, pimg->timetag,-1);
      pchar = g_utf8_collate_key (pimg->time,-1);
      gtk_tree_store_set(treestore,&iter2,STRINGUTF8_COL, pchar,-1);      
      gtk_tree_store_set(treestore,&iter2,ISFILENAME_COL, FALSE,-1);
      g_free(pchar);
     }else{ 
      //look fou timetag of the input image
      if(bout){
       pchar = strstr(pimg->fname,".out.tif");
       if(!pchar) strstr(pimg->fname,".out.TIF");
       g_string_insert_len(saux,0,pimg->fname,pchar-(pimg->fname));
       g_string_truncate(saux,pchar-(pimg->fname));
       ttag = GetIMGTimeTag(bf,saux->str);
       if(ttag<0) ttag = pimg->timetag;
      }else
       ttag = pimg->timetag;
      
      GetClosestTimeTagTreeIter(ttag, &iter1, &iter2);
    } 
   }
   

  gtk_tree_store_append(treestore, &iter3, &iter2);
  gtk_tree_store_set(treestore,&iter3,STRING_COL , pimg->fname,-1);
  gtk_tree_store_set(treestore,&iter3,TIMETAG_COL, pimg->timetag,-1);
  gtk_tree_store_set(treestore,&iter3,ISFILENAME_COL, TRUE,-1);
  pchar2=g_utf8_casefold(pimg->fname,-1);
  pchar = g_utf8_collate_key_for_filename (pchar2,-1);
  gtk_tree_store_set(treestore,&iter3,STRINGUTF8_COL, pchar,-1);
  g_free(pchar);
  g_free(pchar2);
  
  n = n->next;
 }
 g_string_free(saux,TRUE);
}




gint GetClosestTimeTagTreeIter(gdouble dtime,
                                GtkTreeIter *parent,
                                GtkTreeIter *iter){
    gboolean valid;
    gint iRes;
    double daux;
    GtkTreeModel *treemodel;
    gint i,icmp,icmp_old;
    gdouble ttag_old;
    
   
   
    treemodel = GTK_TREE_MODEL(treestore); 
    valid = gtk_tree_model_iter_children(treemodel, iter, parent);
    if(valid==0)
     return(0);  //no child nodes yet...return code to create one 
     
    i=0;
    icmp_old=-1;
    ttag_old=-99999.0;
    do{
       gtk_tree_model_get(treemodel, iter,TIMETAG_COL, &daux,-1);
       icmp=(daux>dtime)?1:-1;
       if(dtime==daux || icmp==1){  //exact match or dtime
        break;                      //earlier than TIMETAG;
       }
       i++;
       icmp_old = icmp; 
       ttag_old = daux;
     }while(gtk_tree_model_iter_next (treemodel,iter));
     if( (i<gtk_tree_model_iter_n_children(treemodel,parent)) &&    
         (dtime-ttag_old)<(daux-dtime)){ //previous node is closest
           gtk_tree_model_iter_nth_child(treemodel,iter,parent,i-1);
           gtk_tree_model_get(treemodel, iter,TIMETAG_COL, &daux,-1);
     }
     
     return(iRes);
}


//la funcion devuelve 1 si se encontro el time y se
//fija iter al pos encontrado
//se devuelve 0 si no, en este caso iter es fijado al pos (BF) mas
//cercano ya existente.
//iRes: 1 - exact matchgtk_tree_path_new_from_string
//iRes: 0 - closest match
//iRes:-1 - ohoh.
gint GetClosestTimeTreeIter(gchar *stime,
                                GtkTreeIter *parent,
                                GtkTreeIter *iter){
    gboolean valid;
    gint iRes = 0;
    gchar *pchar1;
    GtkTreeModel *treemodel;
    gint i,icmp,icmp_old;
    gchar *saux;
    
   

    treemodel = GTK_TREE_MODEL(treestore); 
    valid = gtk_tree_model_iter_children(treemodel, iter, parent);
    if(valid==0)
     return(0);  //no child nodes yet...return code to create one 
    
    saux = g_utf8_collate_key(stime,-1); 
    i=0;
    icmp_old=0;
    do{
       gtk_tree_model_get(treemodel, iter,STRINGUTF8_COL, &pchar1,-1);
       icmp=strcmp(pchar1,saux);
       g_free(pchar1);
       if(icmp==0){  //exact match
        iRes=1;
        break;
       }
       if(icmp>0){
        if(i==0){ //problems: stime < first node
         iRes=-1;
         break;
        }else{
          iRes=0;
          gtk_tree_model_iter_nth_child(treemodel, iter, parent, i-1);
          break; 
        } 
       }
       i++;
       icmp_old = icmp; 
     }while(gtk_tree_model_iter_next (treemodel,iter));
     g_free(saux);
     return(iRes);
}

gboolean GetTreeIter(gchar *saux, GtkTreeIter *iter, GtkTreeIter *parent){
    gboolean valid,bFound=FALSE;
    gchar *pchar1;
    GtkTreeModel *treemodel;
     
    treemodel = GTK_TREE_MODEL(treestore); 
    valid = gtk_tree_model_iter_children(treemodel, iter, parent);
    while(valid){
      gtk_tree_model_get(treemodel, iter,STRING_COL, &pchar1,-1);
      if(strcmp(pchar1,saux)==0){
        bFound = TRUE;
        break;
      }  
      valid = gtk_tree_model_iter_next (treemodel,iter);
      g_free(pchar1);
    }
    return(bFound);
}

// Returns a list of  (IMG *),
// matching *%s[%s.]*tif or *%s[%s.]*TIF patterns
GList *GetFileNames(gchar *ppath, gchar *pattern){
 IMG *pimg;
 GtkWidget   *w;
 GDir        *gdir;
 GString     *saux = g_string_new(NULL);
 const gchar *pfname;
 gchar       *pchar1=NULL,
             *pchar2=NULL,
             *pchar3=NULL;
 GError      *error = NULL;
 GList       *pList = NULL;
 GPatternSpec *pspec1;
 gboolean    b;
 gint        i;
 GError *err;
 gsize  bwritten;
 gchar  *patUTF8;

 gdir = g_dir_open(ppath,0,&error);
 if (gdir==NULL){
   g_string_printf(saux,"g_dir_open(%s) failed - %s\n",ppath ,error->message);
   message_dialog(GTK_MESSAGE_WARNING,GTK_WINDOW(app1),
                  "Error",saux->str);
   g_error_free(error);
   g_string_free(saux,TRUE);
   return(NULL);
 }else{
  //pspec1 = g_pattern_spec_new(pattern);
  
  patUTF8 = g_locale_to_utf8(pattern,-1,NULL,&bwritten,&err);
  pspec1  = g_pattern_spec_new(patUTF8);
  g_free(patUTF8);  

  while((pfname=g_dir_read_name(gdir))) {
    //chequear .out.  
//    if(((g_strrstr(pfname,".out.") ||
//         g_strrstr(pfname,p.output_basename->str))))
//      continue;
	
    //check token file: 
    //looking for "token"
    //i=0; while(pfname[i++]!='\0');
    //b = g_pattern_match(pspec1,i,pfname,NULL);  //ojo UTF-8 Multiplataforma?!
    patUTF8 = g_locale_to_utf8(pfname,-1,NULL,&bwritten,&err);
    i=0; while(patUTF8[i++]!='\0');
    b = g_pattern_match(pspec1,i,patUTF8,NULL);  //ojo UTF-8 Multiplataforma?!
    g_free(patUTF8);

    if(b){
     pchar1 = NULL;
     pchar2 = NULL;
     pimg = (IMG *)malloc(sizeof(IMG));
     pimg->fname=NULL;
     pimg->pos  =NULL;
     pimg->time=NULL;     
     pimg->fname_utf=NULL;
     pimg->pos_utf  =NULL;
     pimg->time_utf =NULL;
     pimg->timetag = -1;
     
     pimg->fname = (gchar *) malloc(sizeof(gchar)*i);
     strcpy(pimg->fname,pfname);
     pimg->fname_utf = g_utf8_collate_key_for_filename(pimg->fname,-1);
     
     if(p.bflat)
      if(strstr(pfname,p.flat_token->str)){
       goto jump1;
      }
     if(p.bdark)
      if(strstr(pfname,p.dark_token->str)){
       goto jump1;
      }
      
     //hay pos y time?
     w = lookup_widget(inputDialog,"cbPosToken");
     if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))){
      pchar1=g_strrstr(pfname,p.pos_token->str);
     }
     w = lookup_widget(inputDialog,"cbTimeToken");
     if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))){
      pchar2=g_strrstr(pfname,p.time_token->str);
     }
     // pos and not time
     if((pchar1!=NULL) && (pchar2==NULL)){
      pchar2=strstr(pchar1,".");
      pimg->pos     = g_strndup(pchar1,pchar2-pchar1);
      pimg->pos_utf = g_utf8_collate_key_for_filename(pimg->pos,-1);
     }else{
       //no pos and time
       if((pchar1==NULL) && (pchar2!=NULL)){
        pchar1=strstr(pchar1,".");
        pimg->time     = g_strndup(pchar2,pchar1-pchar2);
        pimg->time_utf = g_utf8_collate_key_for_filename(pimg->time,-1);
       }else{
        //pos and time
        if((pchar1!=NULL) && (pchar2!=NULL)){
          pchar3=strstr(pchar1,p.separator->str);
          pimg->pos = g_strndup(pchar1,pchar3-pchar1);
          pimg->pos_utf = g_utf8_collate_key_for_filename(pimg->pos,-1);
          pchar3=strstr(pchar2,".");
          pimg->time = g_strndup(pchar2,pchar3-pchar2);
          pimg->time_utf = g_utf8_collate_key_for_filename(pimg->time,-1);
        }
       }
      }
      jump1:
      pList = g_list_append(pList,pimg); 
     } 
    } 
  
  g_dir_close(gdir); 
  g_pattern_spec_free(pspec1);
 }
 return(pList);
}


/* This is not pretty. Of course you can also use a
*  separate compare function for each sort ID value */
gint  sort_func_filename (GtkTreeModel *model,
                          GtkTreeIter  *a,
                          GtkTreeIter  *b,
                          gpointer      userdata){
    gint ret = 0;
    gchar *name1, *name2;
    gboolean b1,b2;
    
    gtk_tree_model_get(model, a, ISFILENAME_COL, &b1, -1);
    gtk_tree_model_get(model, b, STRINGUTF8_COL, &b2, -1);
    
    if(b1 || b2){
     return(-1);  
    }   
      
    return(strcmp(name1,name2));  
      
  }
  


//////////////////////////// ImageView Functions //////////////////
void LoadImages (gchar *fname)
{ 
    GdkPixbufAnimation *anim;
    gboolean is_image;
    GString *path = g_string_new(NULL);
    gchar *pchar=NULL, *pchar1=NULL; 
    GtkWidget *w;
    gchar *fl_path  =p.fl_path->str, 
          *flat_path=p.flat_path->str, 
          *dark_path=p.dark_path->str;
    
    
    if((pchar=g_strrstr(fname,"TEST_"))){
     pchar = p.tmp_path->str;
     g_string_assign(path,pchar);
     fname = &(fname[5]);
    }
    else{
     if(p.same_path){
      //fl_path = p.bf_path->str;
      //flat_path = p.bf_path->str;
      //dark_path = p.bf_path->str;
      g_string_assign(path,p.bf_path->str);
     }else{
      //build absolute path from name
      if((pchar=strstr(fname,p.flat_token->str)) && p.bflat){
         g_string_assign(path,flat_path);
      }else{ 
       if((pchar=strstr(fname,p.dark_token->str)) && p.bdark){
         g_string_assign(path,dark_path);
       }else{ 
        if((pchar=strstr(fname,p.bf_token->str)) &&
           (strcmp(p.bf_token->str,"")!=0)){
         if(pchar==fname){ //make sure the match is at the beginning
          g_string_assign(path,p.bf_path->str);
         } 
        }else{
         if((pchar=strstr(fname,p.fl_token->str)) && 
            (strcmp(p.fl_token->str,"")!=0) ){
          if(pchar==&(fname[1])){ //make sure the match is at the beginning
           g_string_assign(path,fl_path);
          } 
         }
        }  
       }
      }  
     }
    } 
    pchar = g_build_filename(path->str,fname,NULL); 
    g_string_assign(path,pchar);
    g_free(pchar);
    
    
    anim = gdk_pixbuf_animation_new_from_file (path->str,NULL);
    if (!anim)
    {
        printf ("No anim!\n");
        return;
    }
    gtk_anim_view_set_anim (view, anim);
    g_object_unref (anim);

    char *basename = g_path_get_basename (path->str);
    gtk_window_set_title (GTK_WINDOW(app1), basename);
 //   push_image_info (basename, anim);
    g_free (basename);
    g_string_free(path,TRUE);
    //gtk_action_group_set_sensitive (image_group, TRUE);
//
//    /* Only active the transform_group if the loaded object is a single
//       image -- transformations cannot be applied to animations. */
   is_image = gdk_pixbuf_animation_is_static_image (anim);
//    gtk_action_group_set_sensitive (transform_group, is_image);
  
   w = lookup_widget(app1,"tbZoomIn");
   gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);
   w = lookup_widget(app1,"tbZoomOut");
   gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);
   w = lookup_widget(app1,"tbZoom1");
   gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);
   w = lookup_widget(app1,"tbZoomFit");
   gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);
   
}


gint GetTreeDepth(GtkTreeView *tv){
 GtkTreeIter iter,parent;
 GtkTreeIter *piter,*pparent;
 gint i=0;
 gboolean valid;
 
 piter   = &iter;
 pparent = &parent;
 valid = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(treestore),pparent);
 while(valid){
  i++;
  valid= gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(treestore),piter,pparent,0); 
  pparent = piter;
 }
 return(i);
}

///////////////////////////// Popup Menu //////////////////////////
void  view_popup_menu_onDoTest (GtkWidget *menuitem, gpointer userdata){
  GtkWidget *w;
  
  gboolean b;
   //Setup de segmentDialog in its 'Test' form. 
  w = lookup_widget(segmentDialog,"bRun");
  g_object_get(G_OBJECT(w),"visible",&b,NULL);
  if(b==TRUE){
   gtk_widget_hide(GTK_WIDGET(segmentDialog));
   g_object_set(G_OBJECT(w),"visible",FALSE,NULL);
//   w = lookup_widget(segmentDialog,"frameAlignement");
//   g_object_set(G_OBJECT(w),"visible",FALSE,NULL);  
   w = lookup_widget(segmentDialog,"bTest");
   g_object_set(G_OBJECT(w),"visible",TRUE,NULL);
  } 
  gtk_widget_show(GTK_WIDGET(segmentDialog));
  
  //Bug: donot know why but if not set explicitely it does not appear
//   w = lookup_widget(segmentDialog,"frameAlignement");
//   g_object_set(G_OBJECT(w),"visible",TRUE,NULL);
}


void  view_popup_menu_onDoRun (GtkWidget *menuitem, gpointer userdata)
{
  GtkWidget *w;
  gboolean b;
  
  //Setup de segmentDialog in its 'Run' form.
  w = lookup_widget(segmentDialog,"bRun");
  g_object_get(G_OBJECT(w),"visible",&b,NULL);
  if(b==FALSE){
   gtk_widget_hide(GTK_WIDGET(segmentDialog));
   g_object_set(G_OBJECT(w),"visible",TRUE,NULL);
//    w = lookup_widget(segmentDialog,"frameAlignement");
//    g_object_set(G_OBJECT(w),"visible",TRUE,NULL);   
   w = lookup_widget(segmentDialog,"bTest");
   g_object_set(G_OBJECT(w),"visible",FALSE,NULL);
  } 
  gtk_widget_show(GTK_WIDGET(segmentDialog));
  
  //Bug: donot why but if not set explicitely it does not appear  
//   w = lookup_widget(segmentDialog,"frameAlignement");
// //   g_object_set(G_OBJECT(w),"visible",TRUE,NULL);
}


void  view_popup_menu (GtkWidget *treeview, GdkEventButton *event, gpointer userdata)
  {
    GtkWidget *menu, *menuitem;
    GtkWidget *w;
    gboolean b1,b2;
    GtkTreeSelection *selection;
    GtkTreeView *tv;
    GtkTreeModel *treemodel;
    GtkTreeIter  iter, child, *piter;
    gchar *fname;
    gint breturn;
    
    //22/10/2008: check if launching from root child
    //do not show popup menu on dark or flat entries.
    w = lookup_widget(app1,"treeview1");
    tv = GTK_TREE_VIEW(w);
    treemodel = gtk_tree_view_get_model(tv);
    selection  = gtk_tree_view_get_selection(GTK_TREE_VIEW(tv));
    if(!gtk_tree_selection_get_selected(selection,&treemodel,&iter))
     return;
    gtk_tree_model_get(treemodel,&iter,0,&fname,-1);
    breturn  = FALSE;
    if(p.bdark)
     if(strstr(fname,p.dark_token->str) ||
        strstr(fname,"Backgd Illumination")){
         breturn = TRUE;
     }
    if(p.bflat)
     if(strstr(fname,p.flat_token->str) ||
        strstr(fname,"Uneven Illumination")){
         breturn = TRUE;
     }
    g_free(fname); 
    if(breturn){
     return;
    } 
     
    
    menu = gtk_menu_new();

    w  = lookup_widget(segmentDialog,"bTest");
    g_object_get(GTK_WIDGET(segmentDialog),"visible",&b1,NULL);
    g_object_get(GTK_WIDGET(w),"visible",&b2,NULL);
    if(!b1 || (b1 && !b2)){
     menuitem = gtk_menu_item_new_with_label("Show Test Dialog");
     g_signal_connect(menuitem, "activate",
                     (GCallback) view_popup_menu_onDoTest, treeview);
     gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
    }
    if(!b1 || (b1 && b2)){
     menuitem = gtk_menu_item_new_with_label("Show Run CellID Dialog");
     g_signal_connect(menuitem, "activate",
                     (GCallback) view_popup_menu_onDoRun, treeview);
     gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
    }
    gtk_widget_show_all(menu);

    /* Note: event can be NULL here when called from view_onPopupMenu;
     *  gdk_event_get_time() accepts a NULL argument */
    gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL,
                   (event != NULL) ? event->button : 0,
                   gdk_event_get_time((GdkEvent*)event));
  }


gboolean
  view_onPopupMenu (GtkWidget *treeview, gpointer userdata)
  {
    view_popup_menu(treeview, NULL, userdata);

    return TRUE; /* we handled this */
  }


////////////////////  Launch CellID  ////////////////////////////

void TestCellID(){ 
    GtkTreeView *treeview ;
    GtkTreeSelection *selection;
    GtkTreeIter    iter,  iterAux, iterAux2;
    GtkTreeIter  *piter, *piterAux;
    gint idepth, idepthIter,i;
    gboolean valid;
    char *pchar1, *basename,*pchar2;
    GtkTreeModel *treemodel;
    GtkWidget *w;
    GtkTreePath *gpath;
    GString *saux=g_string_new(NULL);
    
    //Check existence and availability of cellid exec  
    if(!CheckExecFile()){
     gtk_widget_show(GTK_WIDGET(setupDialog));
     return;
    }
    
    w = lookup_widget(app1,"treeview1");
    treeview = GTK_TREE_VIEW(w);
    
    treemodel = gtk_tree_view_get_model(treeview);
    
    selection  = gtk_tree_view_get_selection(GTK_TREE_VIEW(treeview));
    if(!gtk_tree_selection_get_selected(selection,&treemodel,&iter))
     return;

    idepth = 1;
    if(p.bpos ) idepth++;
    if(p.btime) idepth++;
    idepthIter = gtk_tree_store_iter_depth(treestore,&iter);
    if(idepthIter>idepth){
      gtk_tree_model_iter_parent(GTK_TREE_MODEL(treestore),&iterAux,&iter);
      iter = iterAux;
      idepthIter--;
    }
    if(idepthIter!=idepth){
     message_dialog(GTK_MESSAGE_WARNING,
                     GTK_WINDOW(app1),
                     "Warning",
                     "Select a bright field image to proceed");

     return;
    } 
    gtk_tree_model_iter_parent(GTK_TREE_MODEL(treestore),&iterAux,&iter);

    #ifdef MYDEBUG
    gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iterAux,0,&pchar1,-1);
    g_free(pchar1);
    #endif
    
    gtk_widget_show_now(GTK_WIDGET(calculatingDialog));
    
    //valid = LaunchCellID(piter,TRUE); //testmode=TRUE
    valid = LaunchCellID(&iterAux,TRUE,SET_RUN); //testmode=TRUE
    
    gtk_widget_hide(GTK_WIDGET(calculatingDialog));
    
    if(!valid) return;
    
    //insert and select testnode
    //if(GetBFTreeIter(&iter, piter)){
    if(GetBFTreeIter(&iter, &iterAux)){
       gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iter,0,&basename,-1);
       //g_string_printf(saux,"%sTEST%s.out.tif",p.bf_token->str,&(basename[2]));
       g_string_printf(saux,"TEST_%s.out.tif",basename);
       g_free(basename);
       if(gtk_tree_model_iter_has_child(GTK_TREE_MODEL(treestore),&iter))
         gtk_tree_model_iter_children(GTK_TREE_MODEL(treestore),&iterAux,&iter);
       else  
         gtk_tree_store_append(treestore,&iterAux,&iter);
       gtk_tree_store_set(treestore, &iterAux,STRING_COL,saux->str ,-1);
       gtk_tree_store_set(treestore, &iterAux,ISFILENAME_COL,FALSE ,-1);
       
       
       gpath = gtk_tree_model_get_path(GTK_TREE_MODEL(treestore),&iterAux);
       w = lookup_widget(app1,"treeview1");
       gtk_tree_view_expand_to_path(GTK_TREE_VIEW(w),gpath);
       gtk_tree_selection_select_iter(GTK_TREE_SELECTION(selection),&iterAux);
       gtk_tree_path_free(gpath);
       LoadImages(saux->str);       
    }   
    /*
    pchar1 = InsertTestBFNode(piter);
    if(pchar1){
      gpath = gtk_tree_path_new_from_string(pchar1);
      gtk_tree_selection_select_path(GTK_TREE_SELECTION(selection),gpath);
      gtk_tree_path_free(gpath);
    }  */
    g_string_free(saux,TRUE);
}


void RunCellID(){ 
    GtkTreeView *treeview ;
    GtkTreeSelection *selection;
    GtkTreeIter    iter,  iterAux, iterAux2;
    GtkTreeIter  *piter, *piterAux;
    gint idepth, idepthIter,i;
    gboolean valid;
    char *pchar1, *basename,*pchar2;
    GtkTreeModel *treemodel;
    GtkWidget *w;
    GtkTreePath *gpath;
    GString *saux=g_string_new(NULL);
    gint i_type;
    
    //Check existence and availability of cellid exec  
    if(!CheckExecFile()){
     gtk_widget_show(GTK_WIDGET(setupDialog));
     return;
    }
    
    w = lookup_widget(app1,"treeview1");
    treeview = GTK_TREE_VIEW(w);
    
    treemodel = gtk_tree_view_get_model(treeview);
    
    selection  = gtk_tree_view_get_selection(GTK_TREE_VIEW(treeview));
    if(!gtk_tree_selection_get_selected(selection,&treemodel,&iter))
     return;

    idepth = 1;
    if(p.bpos ) idepth++;
    if(p.btime) idepth++;
    idepthIter = gtk_tree_store_iter_depth(treestore,&iter);
    if(idepthIter>idepth){
      gtk_tree_model_iter_parent(GTK_TREE_MODEL(treestore),&iterAux,&iter);
      iter = iterAux;
      idepthIter--;
    }
    
    //SET_RUN.no loop
    //POS_RUN loop over a 1 level loop to generate bf and fl input files 
    //MULTIPOS_RUN iterative POS_RUN 
    if(idepth==3)
     switch(idepthIter){
      case 0:
       iterAux = iter;
       i_type = MULTIPOS_RUN;
       break;
       
      case 1:
       iterAux= iter;
       i_type = POS_RUN;
       break;
       
      case 2: 
        iterAux = iter;
        i_type = SET_RUN; 
        break;

      case 3:
        gtk_tree_model_iter_parent(GTK_TREE_MODEL(treestore),&iterAux,&iter);
        i_type = SET_RUN; 
        break;
     }  
    
    if(idepth==2)
     switch(idepthIter){
      case 0:
       iterAux = iter;
       i_type = MULTIPOS_RUN;
       break;
       
      case 1: 
        iterAux = iter;
        i_type = SET_RUN; 
        break;
       
      case 2:
        gtk_tree_model_iter_parent(GTK_TREE_MODEL(treestore),&iterAux,&iter);
        i_type = SET_RUN; 
        break;
     }  

    #ifdef MYDEBUG
    gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iterAux,0,&pchar1,-1);
    g_free(pchar1);
    #endif
    
    gtk_widget_show_now(GTK_WIDGET(calculatingDialog));
    
    b_coming_from_multipos = FALSE;
    if(i_type==MULTIPOS_RUN){
     b_coming_from_multipos = TRUE;
     w = lookup_widget(calculatingDialog,"label91");
     iter = iterAux;
     valid = gtk_tree_model_iter_children(GTK_TREE_MODEL(treestore),&iterAux, &iter); 
     while(valid){
       gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iterAux,0,&pchar1,-1);
       gtk_label_set_label(GTK_LABEL(w),pchar1);
       g_free(pchar1);
       //if two factors (pos + time) for each pos run multitime...POS_RUN
       if(idepth==3)
        LaunchCellID(&iterAux,FALSE,POS_RUN);
       //if one factor (pos) run each pos indendently 
       if(idepth==2) 
        LaunchCellID(&iterAux,FALSE,SET_RUN);
       valid=gtk_tree_model_iter_next(GTK_TREE_MODEL(treestore),&iterAux);
     }
    }else{
     valid = LaunchCellID(&iterAux,FALSE,i_type); //testmode=FALSE
    } 
    
    gtk_widget_hide(GTK_WIDGET(calculatingDialog));

    LoadImageNames();
    g_string_free(saux,TRUE);
}


//insert test output as a child
//returns the string path of the testRow
// or NULL if something bad happened
gchar *InsertTestBFNode__A(GtkTreeIter *piter){
    GString *saux = g_string_new(NULL);
    GPatternSpec *pspec1;
    gboolean valid;
    GtkTreeIter iterAux,iterAux2;
    gint     icount=0,i;
    gint     iBF  =-1;
    gboolean bTest=FALSE, b;
    gchar *basename,*pchar;
        
        
    g_string_printf(saux,"%s*",p.bf_token->str);
    pspec1 = g_pattern_spec_new(saux->str);
    icount = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(treestore),piter);
    for(i=0; i<icount;i++){
      gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(treestore),&iterAux,piter,i);
      gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iterAux,0,&basename,-1);
      if(!(g_strrstr(basename,".out.") && g_strrstr(basename,p.output_basename->str))){
       b = g_pattern_match_string(pspec1,basename);  //ojo UTF-8 Multiplataforma?!
       if(b){ 
        iBF = i ;
        //check if test node already attached to BF node
        if(gtk_tree_model_iter_has_child(GTK_TREE_MODEL(treestore),&iterAux))
         bTest=TRUE;
        break;
       } 
      }
      g_free(basename); 
    }
     
     
    if(iBF>-1){
     if(!bTest){
      gtk_tree_store_append(treestore, &iterAux2, &iterAux);
      g_string_printf(saux,"%s.%s.tif",basename,p.output_basename->str);
      gtk_tree_store_set(treestore   , &iterAux2,STRING_COL,saux->str ,-1);
      gtk_tree_store_set(treestore   , &iterAux2,ISFILENAME_COL, FALSE,-1);
     }else
      gtk_tree_model_iter_children(GTK_TREE_MODEL(treestore),&iterAux2,&iterAux);
      
     pchar=gtk_tree_model_get_string_from_iter(GTK_TREE_MODEL(treestore),&iterAux2);
    }else{
     message_dialog(GTK_MESSAGE_ERROR,
                     GTK_WINDOW(app1),
                     "Error",
                     "No BF found to attach test image (weird)");
     pchar=NULL;
    }
   
   g_pattern_spec_free(pspec1);
   g_string_free(saux,TRUE); 
   return(pchar);
}


gchar *InsertTestBFNode(GtkTreeIter *piter){
    GString *saux = g_string_new(NULL);
    GPatternSpec *pspec1;
    gboolean valid;
    GtkTreeIter iterAux,iterAux2;
    gint     icount=0,i;
    gint     iBF  =-1;
    gboolean bTest=FALSE, b;
    gchar *basename,*pchar;
        
        
    g_string_printf(saux,"%s*",p.bf_token->str);
    pspec1 = g_pattern_spec_new(saux->str);
    
    icount = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(treestore),piter);
    for(i=0; i<icount;i++){
      gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(treestore),&iterAux,piter,i);
      gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iterAux,0,&basename,-1);
      if(!(g_strrstr(basename,".out.") && g_strrstr(basename,p.output_basename->str))){
       b = g_pattern_match_string(pspec1,basename);  //ojo UTF-8 Multiplataforma?!
       if(b){ 
        iBF = i ;
        //check if test node already attached to BF node
        if(gtk_tree_model_iter_has_child(GTK_TREE_MODEL(treestore),&iterAux))
         bTest=TRUE;
        break;
       } 
      }
      g_free(basename); 
    }
     
     
    if(iBF>-1){
     if(!bTest){
      gtk_tree_store_append(treestore, &iterAux2, &iterAux);
      g_string_printf(saux,"%s.%s.tif",basename,p.output_basename->str);
      gtk_tree_store_set(treestore   , &iterAux2,STRING_COL,saux->str ,-1);
      gtk_tree_store_set(treestore   , &iterAux2,ISFILENAME_COL, FALSE,-1);
     }else
      gtk_tree_model_iter_children(GTK_TREE_MODEL(treestore),&iterAux2,&iterAux);
      
     pchar=gtk_tree_model_get_string_from_iter(GTK_TREE_MODEL(treestore),&iterAux2);
    }else{
     message_dialog(GTK_MESSAGE_ERROR,
                     GTK_WINDOW(app1),
                     "Error",
                     "No BF found to attach test image (weird)");
     pchar=NULL;
    }
   
   g_pattern_spec_free(pspec1);
   g_string_free(saux,TRUE); 
   return(pchar);
}


void  WriteTXTfiles(GtkTreeIter *parent,gboolean btest,gint itype){
 GString *saux = g_string_new(NULL);
 GtkTreeIter iter, iter0;
 gint i=0;
 gboolean valid, valid1;
 
 switch(itype){
  case SET_RUN:
    WriteTXTfiles_Set(parent,"w",btest);
    break;
    
  case POS_RUN: 
  case MULTIPOS_RUN: 
    valid=gtk_tree_model_iter_children(GTK_TREE_MODEL(treestore),&iter, parent);
    i=0;
    while(valid){
     if(i==0)
      WriteTXTfiles_Set(&iter,"w",FALSE);
     else 
      WriteTXTfiles_Set(&iter,"a",FALSE);   
     i++;
     valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(treestore),&iter) ;
    }
    break;
 }
}


gboolean LaunchCellID(GtkTreeIter *pIterSet,
                          gboolean btest,
                          gint itype){	
  gchar *argv[5];
  GString *saux0=g_string_new(NULL),
          *saux1=g_string_new(NULL),
          *saux2=g_string_new(NULL),
          *saux3=g_string_new(NULL),
          *saux4=g_string_new(NULL),
          *saux5=g_string_new(NULL);  
  GString *saux=g_string_new(NULL);
  gchar  *sstderr, *sstdout, *fname = NULL, *fname_new = NULL;
  GError *err = NULL;
  gint exit_status;
  gboolean b,valid;
  GtkTextBuffer *tbuf;
  GtkTextIter   start, end;
  GtkWidget *w;
  GtkTreeIter iter;
  gchar *basename,*pchar1,*pchar2,*pcharTmpDir=NULL;
  gchar *contents;
  gsize length;
  gint iaux;

  if(btest){ //copy bf image to tmp file
   valid = TRUE;
   b=GetBFTreeIter(&iter, pIterSet);
   if(b){
    pcharTmpDir = p.tmp_path->str;
    gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iter,0,&basename,-1);
    fname     = g_build_filename(p.bf_path->str,basename,NULL);
    fname_new = g_build_filename(pcharTmpDir,basename,NULL);
    b = g_file_get_contents(fname,  &contents, &length, &err);
    if(!b){
     g_string_printf(saux,"Error reading %s:\nEcode:%d\nEcode:%s\n",
                          fname,err->code,err->message);
     valid=FALSE;
    }
    b = g_file_set_contents(fname_new, contents, length, &err);
    if(!b){
     g_string_printf(saux,"Error reading %s:\nEcode:%d\nEcode:%s\n",
                          fname_new,err->code,err->message);
     valid=FALSE;
    }
    g_free(fname_new);
   }else{
    valid=FALSE;
    g_string_assign(saux,"Problems getting bright field info");
   }
   if(!valid){
    message_dialog(GTK_MESSAGE_ERROR,
                   GTK_WINDOW(app1),
                   "Error",
                    saux->str);
    return(FALSE);                
   }                  
  }
  
  //set p.output_path to bf_path/pos_token dir
  //create dir if necessary
  if(!btest){
    gtk_tree_model_get(GTK_TREE_MODEL(treestore),pIterSet,0,&basename,-1); 
    //if paw output is selected and we are performing a root level analysis
    //then export results into posX folders nomatter postoken
    if(p.bPAW && b_coming_from_multipos){
     pchar1 = strstr(basename,p.pos_token->str);
//     if(p.btime)
//       pchar2 = strstr(basename,p.separator->str);
//     else
//       pchar2 = strstr(basename,".");
     pchar1=pchar1+p.pos_token->len;
//     pchar2[0]='\0';
     iaux = atoi(pchar1);
     g_string_printf(saux,"pos%d",iaux);
     fname = g_build_filename(p.bf_path->str,saux->str,NULL); 
    }else{
     fname = g_build_filename(p.bf_path->str,basename,NULL); 
    } 
    if(!g_file_test(fname,G_FILE_TEST_IS_DIR))
      g_mkdir(fname,0755);
    g_string_assign(p.output_path,fname);
    g_free(basename);
    g_free(fname);
  }
  

  WriteTXTfiles(pIterSet,btest,itype);//fl_vcellid.txt , bf_vcellid.txt
  SaveParams(btest);                 //"parameters_vcellid_%s.txt",p.output_basename->str
 
  g_string_assign(saux0,p.exec_file->str);
  #ifdef G_OS_WIN32
   g_string_append(saux0,"'");
   g_string_prepend(saux0,"'");
  #endif
  argv[0] = saux0->str;
  #ifdef G_OS_WIN32
   g_strdelimit(argv[0],"/",G_DIR_SEPARATOR);
   g_strdelimit(saux0->str,"/",G_DIR_SEPARATOR);
  #endif
  
  
  if(btest)
    fname = g_build_filename(pcharTmpDir,"bf_vctest.txt",NULL);
  else  
    fname = g_build_filename(p.output_path->str,"bf_vcellid.txt",NULL);
   g_string_assign(saux,fname);
  #ifdef G_OS_WIN32
   g_string_append(saux,"'");
   g_string_prepend(saux,"'");
  #endif  
  g_string_printf(saux1,"-b %s",saux->str);
  g_free(fname);
  argv[1] = saux1->str;
  
  if(btest)
    fname = g_build_filename(pcharTmpDir,"fl_vctest.txt",NULL);
  else  
    fname = g_build_filename(p.output_path->str,"fl_vcellid.txt",NULL);
   g_string_assign(saux,fname);
  #ifdef G_OS_WIN32
   g_string_append(saux,"'");
   g_string_prepend(saux,"'");
  #endif  
  g_string_printf(saux2,"-f %s",saux->str);



  if(!btest && p.bnucleus){
   fname = g_build_filename(p.output_path->str,"nuc_vcellid.txt",NULL);
   #ifdef G_OS_WIN32
    g_string_printf(saux," -3 '%s'",fname);
   #else
    g_string_printf(saux," -3 %s",fname);
   #endif
   g_string_append(saux2,saux->str);
   argv[2] = saux2->str;
  }else{ 
    argv[2] = saux2->str;
  }
  g_free(fname);
   
  
  
  if(btest){
   g_string_printf(saux,"parameters.txt");
   fname = g_build_filename(pcharTmpDir,saux->str,NULL);
  } 
  else{
   g_string_printf(saux,"parameters_vcellid_%s.txt",p.output_basename->str);
   fname = g_build_filename(p.output_path->str,saux->str,NULL);
  } 
  g_string_assign(saux,fname);
  #ifdef G_OS_WIN32
   g_string_append(saux,"'");
   g_string_prepend(saux,"'");
  #endif  
  g_string_printf(saux4,"-p %s",saux->str);
  g_free(fname);
  argv[3] = saux4->str;
    
  
  if(btest){
   fname = g_build_filename(pcharTmpDir,p.output_basename->str,NULL); 
   g_string_assign(saux,fname);
  #ifdef G_OS_WIN32
   g_string_append(saux,"'");
   g_string_prepend(saux,"'");
  #endif  
   g_string_printf(saux3,"-o %s",saux->str);
  }else{
    fname = g_build_filename(p.output_path->str,p.output_basename->str,NULL); 
    g_string_assign(saux,fname);
   #ifdef G_OS_WIN32
    g_string_append(saux,"'");
    g_string_prepend(saux,"'");
   #endif  
    g_string_printf(saux3,"-o %s",saux->str);
  } 
  argv[4] = saux3->str;
  g_free(fname); fname = NULL;
    
  if(p.bflat){  
   if(btest)
    fname = g_build_filename(pcharTmpDir,"flat_vctest.txt",NULL);
   else  
    fname = g_build_filename(p.output_path->str,"flat_vcellid.txt",NULL);
   g_string_assign(saux,fname);
   #ifdef G_OS_WIN32
   g_string_append(saux,"'");
   g_string_prepend(saux,"'");
   #endif  
   g_string_printf(saux5," -t %s",saux->str);
   g_string_append(saux2,saux5->str);
   g_free(fname);fname=NULL;
  }

  if(p.bdark){  
   if(btest)
    fname = g_build_filename(pcharTmpDir,"dark_vctest.txt",NULL);
   else  
    fname = g_build_filename(p.output_path->str,"dark_vcellid.txt",NULL);
   g_string_assign(saux,fname);
   #ifdef G_OS_WIN32
   g_string_append(saux,"'");
   g_string_prepend(saux,"'");
   #endif  
   g_string_printf(saux5," -d %s",saux->str);
   g_string_append(saux2,saux5->str);
   g_free(fname);fname=NULL;
  }  
    
    
 if(FALSE){
  b = g_spawn_sync(NULL,//p.output_path->str, 
                   argv, 
                   NULL,
                   //(GSpawnFlags)(G_SPAWN_SEARCH_PATH |G_SPAWN_FILE_AND_ARGV_ZERO),
                   (GSpawnFlags)(G_SPAWN_FILE_AND_ARGV_ZERO),
                   NULL, 
                   NULL, 
                   &sstdout, 
                   &sstderr, 
                   NULL,//&exit_status, 
                   &err)         == FALSE;
  }else{
   g_string_printf(saux,"%s %s %s %s %s",
                   saux0->str,saux1->str,saux2->str,saux3->str,saux4->str);
   if(btest){
     g_string_append(saux," --bf-as-fl");

     //get BF filename and set a new name to rename cellid output file
     /*if(GetBFTreeIter(&iter, pIterSet)){
       gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iter,0,&basename,-1);
       g_string_printf(saux2,"%s.out.tif",basename);
       fname     = g_build_filename(p.bf_path->str,saux2->str,NULL);       
       g_string_printf(saux2,"%sTEST%s.out.tif",p.bf_token->str,&(basename[2]));
       fname_new = g_build_filename(p.bf_path->str,saux2->str,NULL);
       g_free(basename); 
       //if there is an out file, temporarily rename it
       if(g_file_test(fname,G_FILE_TEST_EXISTS))
          g_rename(fname,"tmp.tif");
     } */
   }  
   #ifdef MYDEBUG
    g_printf("%s\n",saux->str);
   #endif
   b = g_spawn_command_line_sync(saux->str,
                                 &sstdout,
                                 &sstderr,
                                 NULL,
                                 &err);
  }
  
  if(!b){
     g_string_printf(saux,"Error code: %d\n%s",err->code,err->message);
     message_dialog(GTK_MESSAGE_ERROR,
                     GTK_WINDOW(app1),
                     "Error",
                     saux->str);
     g_error_free(err);
  }else{  
     g_printf("stdout////////////////\n%s\n",sstdout);
     g_printf("stderr////////////////\n%s\n",sstderr);
     
     //if testing, rename output file
    /* if(btest){
       g_rename(fname,fname_new);
       g_rename("tmp.tif",fname); 
      } 
     */
     
     w = lookup_widget(app1,"txtViewStdout");
     tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(w));
     gtk_text_buffer_get_bounds( tbuf, &start, &end);
     gtk_text_buffer_delete(tbuf, &start, &end);
     gtk_text_buffer_set_text( tbuf, sstdout,-1);

     w = lookup_widget(app1,"txtViewStderr");
     tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(w));
     gtk_text_buffer_set_text( tbuf, sstderr,-1);
     
  }
  


      
  g_string_free(saux,TRUE);
  g_string_free(saux1,TRUE);
  g_string_free(saux2,TRUE);
  g_string_free(saux3,TRUE);
  g_string_free(saux4,TRUE);
  g_string_free(saux5,TRUE);
  g_free(sstdout);
  g_free(sstderr);
  return(b); 
}

// If btest==TRUE run a 'fake' cellid session
// Provide an empty fp_list.txt, and a single entry bf_list.txt
void WriteTXTfiles_Set(GtkTreeIter *pIterSet,gchar *wmode,gboolean btest){ 
  FILE *fpFL=NULL, 
       *fpBF=NULL,
       *fp3 =NULL,
       *fpFlat = NULL,
       *fpDark = NULL;
  gint numFL,num3,i;
  GString *saux=g_string_new(NULL);
  gchar *fname,*basename,*pchar,*pchar2,*fl_path,*fname_nuc=NULL;
  GtkTreeIter iterAux, iter;
  GPatternSpec *pspec1;
  gboolean valid,b,bb,valid2;
  gchar *contents;
  gsize  length;
  GError *err = NULL;
  

  
  
  //first, FL files 
  //write FL input  file
  if(btest){
   pchar = p.tmp_path->str;
   fname = g_build_filename(pchar,"fl_vctest.txt",NULL);
   fpFL = fopen(fname,"w");
   fclose(fpFL);
   numFL = 1;
   num3  = 1;
  }else{
   if(p.same_path)
    fl_path = p.bf_path->str;
   else
    fl_path = p.fl_path->str;
   fname = g_build_filename(p.output_path->str,"fl_vcellid.txt",NULL); 
   fpFL = fopen(fname,wmode);
   if(fpFL==NULL){
     g_string_printf(saux,"%s could not be opened for writing",fname);
     message_dialog(GTK_MESSAGE_ERROR,
                    GTK_WINDOW(app1),
                     "Error",
                     saux->str);
     g_free(fname);
     g_string_free(saux,TRUE);
     return;
   }
   g_free(fname); 
   
   if(p.bnucleus){ 
    fname = g_build_filename(p.output_path->str,"nuc_vcellid.txt",NULL); 
    fp3 = fopen(fname,wmode);
    if(fp3==NULL){
     g_string_printf(saux,"%s could not be opened for writing",fname);
     message_dialog(GTK_MESSAGE_ERROR,
                    GTK_WINDOW(app1),
                     "Error",
                     saux->str);
     g_free(fname);
     g_string_free(saux,TRUE);
     return;
    }
    g_free(fname); 
   }
    
   numFL = 0;num3=0;
   g_string_printf(saux,"?%s*",p.fl_token->str);
   pspec1 = g_pattern_spec_new(saux->str);
   valid = gtk_tree_model_iter_children(GTK_TREE_MODEL(treestore),&iterAux,pIterSet);
   while(valid) {
    gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iterAux,0,&basename,-1);
    valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(treestore),&iterAux);
    //chequear .out.  
    if(((g_strrstr(basename,".out.") || g_strrstr(basename,p.output_basename->str))))
      continue;
    		
    //looking for fluorescence token
    b = g_pattern_match_string(pspec1,basename);  //ojo UTF-8 Multiplataforma?!
    if(b){
     //dont mess with images used to tag the nucleus    		
     if(basename[0]=='3') continue;

     //record florescence filename
     fname = g_build_filename(fl_path,basename,NULL);
     fprintf(fpFL,"%s\n",fname);
     numFL++;

     if(p.bnucleus){
      g_string_insert_len(saux,0,basename,3);
      g_string_set_size(saux,3);

      //we assume first 3 letters are ?FP AND
      //that there is as many fluorescence nuclear tagged images
      //as any other channel (not more, nor less)
      //if fluo file is the selected nuclear channel...proceed
      g_string_printf(saux,"3%s",&(basename[1]));
      fname_nuc = g_build_filename(fl_path,saux->str,NULL);
      if(strstr(basename,p.nucleus_channel->str)){
        //now dump ?FP into nuc image
        bb = g_file_get_contents(fname,&contents, &length, &err);
        if(!bb){
         g_string_printf(saux,"Error reading %s:\nEcode:%d\nEcode:%s\n",
                          fname,err->code,err->message);
         valid=FALSE;
        }
        bb = g_file_set_contents(fname_nuc, contents, length, &err);
        if(!bb){
         g_string_printf(saux,"Error reading %s:\nEcode:%d\nEcode:%s\n",
                          fname_nuc,err->code,err->message);
         valid=FALSE;
        }
      } 
      //add tmp nuclear file into third file list.
      fprintf(fp3,"%s\n",fname_nuc);
      g_free(fname_nuc); 
      num3++;
     }
     g_free(fname);
    }
    g_free(basename);
   }

   fclose(fpFL);
   if(p.bnucleus) fclose(fp3);
   g_pattern_spec_free(pspec1);
  }
  
  //check if every fluorescence file has a nuclear tagged image associated to
  if(p.bnucleus){
   gfloat faux = (numFL % num3);
   if(faux!=0){
     message_dialog(GTK_MESSAGE_WARNING,
                    GTK_WINDOW(app1),
                     "Warning",
                     "Not every fluorescence file has a nuclear tagged\nimage associated to...");
   }
  }
  
  
  //BF
  if(btest)
   fname = g_build_filename(pchar,"bf_vctest.txt",NULL);
  else
   fname = g_build_filename(p.output_path->str,"bf_vcellid.txt",NULL);
  fpBF = fopen(fname,wmode);
  if(fpBF==NULL){
     g_string_printf(saux,"%s could not be opened for writing",fname);
     message_dialog(GTK_MESSAGE_ERROR,
                    GTK_WINDOW(app1),
                     "Error",
                     saux->str);
     g_free(fname);
     g_string_free(saux,TRUE);
     return;
  }

  g_free(fname);
  b=GetBFTreeIter(&iter, pIterSet);
  if(b){
    gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iter,0,&basename,-1);
    if(btest)
     fname = g_build_filename(pchar,basename,NULL);
    else
     fname = g_build_filename(p.bf_path->str,basename,NULL); 
    g_free(basename);
    for(i=0;i<numFL;i++)
      fprintf(fpBF,"%s\n",fname);
    fclose(fpBF);
    g_free(fname);
  }else{

     g_string_printf(saux,"No Bright Field file found!");
     message_dialog(GTK_MESSAGE_ERROR,
                    GTK_WINDOW(app1),
                     "Error",
                     saux->str);
  } 
 
  if(p.bdark){
   if(btest)
    fname = g_build_filename(pchar,"dark_vctest.txt",NULL);
   else 
    fname = g_build_filename(p.output_path->str,"dark_vcellid.txt",NULL);
   fpBF = fopen(fname,"w");
   g_free(fname);
   if(fpBF==NULL){
     g_string_printf(saux,"%s could not be opened for writing",fname);
     message_dialog(GTK_MESSAGE_ERROR,
                    GTK_WINDOW(app1),
                     "Error",
                     saux->str);
     fclose(fpBF);
     g_string_free(saux,TRUE);
     return;
   }
   
   gtk_tree_model_get_iter_first(GTK_TREE_MODEL(treestore),&iterAux);
   valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(treestore),&iterAux);
   while(valid) {
    gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iterAux,0,&fname,-1);
    if(strcmp(fname,"Backgd Illumination")==0){
     g_free(fname);
     valid2 = gtk_tree_model_iter_children(GTK_TREE_MODEL(treestore),&iter,&iterAux);
     while(valid2){
      gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iter,0,&fname,-1);
      if(p.same_path)
       pchar2 = g_build_filename(p.bf_path->str,fname,NULL); 
      else
       pchar2 = g_build_filename(p.dark_path->str,fname,NULL);  
      fprintf(fpBF,"%s\n",pchar2);
      g_free(pchar2);
      g_free(fname);
      valid2 = gtk_tree_model_iter_next(GTK_TREE_MODEL(treestore),&iter);
     } 
     valid = FALSE;
    }else{
     g_free(fname);
     valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(treestore),&iterAux);
    } 
   }
   fclose(fpBF);
  }
  
  if(p.bflat){
   if(btest)
    fname = g_build_filename(pchar,"flat_vctest.txt",NULL);
   else 
    fname = g_build_filename(p.output_path->str,"flat_vcellid.txt",NULL);
   fpBF = fopen(fname,"w");
   g_free(fname);
   if(fpBF==NULL){
     g_string_printf(saux,"%s could not be opened for writing",fname);
     message_dialog(GTK_MESSAGE_ERROR,
                    GTK_WINDOW(app1),
                     "Error",
                     saux->str);
     fclose(fpBF);
     g_string_free(saux,TRUE);
     return;
   }
   
   gtk_tree_model_get_iter_first(GTK_TREE_MODEL(treestore),&iterAux);
   valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(treestore),&iterAux);
   while(valid) {
    gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iterAux,0,&fname,-1);
    if(strcmp(fname,"Uneven Illumination")==0){
     g_free(fname);
     valid2 = gtk_tree_model_iter_children(GTK_TREE_MODEL(treestore),&iter,&iterAux);
     while(valid2){
      gtk_tree_model_get(GTK_TREE_MODEL(treestore),&iter,0,&fname,-1);
      if(p.same_path)
       pchar2 = g_build_filename(p.bf_path->str,fname,NULL); 
      else
       pchar2 = g_build_filename(p.flat_path->str,fname,NULL);  
      fprintf(fpBF,"%s\n",pchar2);
      g_free(pchar2);
      g_free(fname);
      valid2 = gtk_tree_model_iter_next(GTK_TREE_MODEL(treestore),&iter);
     } 
     valid = FALSE;
    }else{
     g_free(fname);
     valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(treestore),&iterAux);
    } 
   }
   fclose(fpBF);
  }
}

// If btest==TRUE run a 'fake' cellid session
gboolean SaveParams(gboolean btest){
 FILE  *fp;
 gchar *fname;
 GString *saux = g_string_new(NULL);
 gchar *pchar;
 
 if(btest){
  pchar = p.tmp_path->str;
  fname = g_build_filename(pchar,"parameters.txt",NULL);
 } 
 else{ 
  g_string_printf(saux,"parameters_vcellid_%s.txt",p.output_basename->str);
  fname = g_build_filename(p.output_path->str,saux->str,NULL);
 } 
 
 fp = fopen(fname,"w");
 if(fp){

   fprintf(fp," max_split_over_minor %.2f\n",p.max_split_over_minor);
   fprintf(fp," max_dist_over_waist %.2f\n",p.max_dist_over_waist);
   if(p.segment_values_t0){
     fprintf(fp," max_split_over_minor_t0 %.2f\n",p.max_split_over_minor);
     fprintf(fp," max_dist_over_waist_t0 %.2f\n",p.max_dist_over_waist);
   }
   
   fprintf(fp," max_pixels_per_cell %d\n",p.max_pixels_per_cell);
   fprintf(fp," min_pixels_per_cell %d\n",p.min_pixels_per_cell);
   
   fprintf(fp," background_reject_factor %.2f\n",p.background_reject_factor);
   fprintf(fp," tracking_comparison %.2f\n",p.tracking_comparison);
   
   if(p.force_nucleus_in_center)
    fprintf(fp," force_nucleus_in_center\n");
    
   if(p.align_cells == CELL_ALIGN_CELL)
    fprintf(fp," align_individual_cells\n");
   if(p.align_cells == CELL_ALIGN_BOUND)
    fprintf(fp," align_individual_cells_boundary\n");
 
  if(p.align_fl == ALIGN_FL_t0)
    fprintf(fp," align_fl_to_first\n");
  if(p.align_fl == ALIGN_FL_BF)
    fprintf(fp," align_fl_to_bf\n");
    
  if(p.image_type==BRIGHTFIELD)
    fprintf(fp," image_type brightfield\n");
  if(p.image_type==CONFOCAL_TRANSMISSION)
    fprintf(fp," image_type confocal_transmisison\n");
  if(p.image_type==METAMORPH_DECONVOLUTION)
    fprintf(fp," image_type metamorph_deconvolution\n");

  //always 'list' for running, but no for testing
  //as we will perform a 'fake' run, with no fluorecent image
  if(btest)
   fprintf(fp," bf_fl_mapping time\n");
  else
   fprintf(fp," bf_fl_mapping list\n");
   
  
  if(p.bPAW)
   fprintf(fp," paw_output\n");
  
  if(p.splitted_fret_image == TRUE){
   fprintf(fp," fret bf_bottom_and_top\n");

   if(p.fret_nucleus == FRET_NUCLEUS_TOP)
    fprintf(fp," fret nuclear_top\n");
   else 
    fprintf(fp," fret nuclear_bottom\n");
  }  
 

 
 
  if(p.treat_brightfield_as_fluorescence)
   fprintf(fp," treat_brightfield_as_fluorescence_also\n");
   
  if(p.bnucleus && !btest)
   fprintf(fp," third_image nuclear_label\n");
   
  fclose(fp);
 }else{
  g_string_printf(saux,"Problems opening %s for saving...",fname);
  message_dialog(GTK_MESSAGE_WARNING,GTK_WINDOW(app1),
                  "Error",saux->str);
  return(FALSE);                
 }
 g_free(fname);
 g_string_free(saux,TRUE);
 return(TRUE);
}


//sets iter to the BF set entry
gboolean GetBFTreeIter(GtkTreeIter *piter, GtkTreeIter *pIterSet){
   GString *saux=g_string_new(NULL);
   GPatternSpec *pspec1;
   gboolean valid,b;
   GtkTreeIter iterAux;
   gchar *basename;
   
   #ifdef MYDEGUB
   gtk_tree_model_get(GTK_TREE_MODEL(treestore),pIterSet,0,&basename,-1);
   g_free(basename);
   #endif
   
   g_string_printf(saux,"%s*",p.bf_token->str);
   pspec1 = g_pattern_spec_new(saux->str);
   valid = gtk_tree_model_iter_children(GTK_TREE_MODEL(treestore),piter,pIterSet);
   while(valid) {
    gtk_tree_model_get(GTK_TREE_MODEL(treestore),piter,0,&basename,-1);
    //chequear .out.  
    if(((g_strrstr(basename,".out.") || g_strrstr(basename,p.output_basename->str)))){
      valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(treestore),piter);
      continue;
    }  
    //check token file: 
    //looking for "token"
    b = g_pattern_match_string(pspec1,basename);  //ojo UTF-8 Multiplataforma?!
    g_free(basename);
    if(b) break;
    valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(treestore),piter);
   }
   g_string_free(saux,TRUE);
   g_pattern_spec_free(pspec1);
   return(valid);
}
