import { BrowserModule } from '@angular/platform-browser';
import { NgModule } from '@angular/core';
import { FormsModule } from '@angular/forms';
import { HttpClientModule } from '@angular/common/http';
import { RouterModule } from '@angular/router';

import { AppComponent } from './app.component';
import { NavMenuComponent } from './nav-menu/nav-menu.component';
import { HomeComponent } from './home/home.component';
import { OptimiseComponent } from './optimise/optimise.component';
import { BarplotComponent } from './barplot/barplot.component';
import { OptimiselossComponent } from './optimiseloss/optimiseloss.component';
import { OptimisegeneralComponent } from './optimisegeneral/optimisegeneral.component';

@NgModule({
  declarations: [
    AppComponent,
    NavMenuComponent,
    HomeComponent,
    OptimiseComponent,
    BarplotComponent,
    OptimiselossComponent,
    OptimisegeneralComponent
  ],
  imports: [
    BrowserModule.withServerTransition({ appId: 'ng-cli-universal' }),
    HttpClientModule,
    FormsModule,
    RouterModule.forRoot([
      { path: '', component: HomeComponent, pathMatch: 'full' },
      {path:'general',component:OptimisegeneralComponent},
      { path: 'etl', component: OptimiseComponent },
      { path: 'loss', component: OptimiselossComponent },
    ])
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
