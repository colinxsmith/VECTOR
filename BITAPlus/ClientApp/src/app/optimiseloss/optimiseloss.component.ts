import { Component, Inject, ElementRef } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import * as d3 from 'd3';

@Component({
  selector: 'app-optimiseloss',
  templateUrl: './optimiseloss.component.html',
  styleUrls: ['./optimiseloss.component.css']
})
export class OptimiselossComponent {
  width = 1000;//width and height for weight graph
  height = 300;
  format = d3.format('0.6f')
  opt: Optimiseloss={}as Optimiseloss;
  constructor(private http: HttpClient, @Inject('BASE_URL') private baseUrl: string, public element: ElementRef) {
    http.get<Optimiseloss>(baseUrl + 'optimise/LOSS').subscribe(result => {
      this.opt=result;
      console.log(this.opt, result);
    }, error => console.error(error));
  }
  over(e: MouseEvent, inout = false) {
    d3.select(e.target as HTMLInputElement & EventTarget).classed('over', inout);
  }
  sendStep() {
    let back = +(d3.select(this.element.nativeElement).select('input.step.g').node() as HTMLInputElement & Event).value;
    console.log(back, this.opt.etLmax);
    this.opt.gamma = back;
    back = +(d3.select(this.element.nativeElement).select('input.step.LU').node() as HTMLInputElement & Event).value;
    this.opt.losSmax = back;
    back = +(d3.select(this.element.nativeElement).select('input.step.LL').node() as HTMLInputElement & Event).value;
    this.opt.losSmin = back;
    this.opt.losSopt = true;
    this.sendData('optimise/LOSS', this.opt)
      .subscribe(ddd => {
        console.log(ddd);
        this.opt = ddd;;
        console.log(this.opt.w);
      }, error => console.error(error));
  }
  sendData(key = 'optimise/LOSS', sendObject: Optimiseloss) {
    const options = {
      headers: new HttpHeaders()
        .set('Content-Type', 'application/json')
    };
    return this.http.post<Optimiseloss>(`${this.baseUrl}${key}`, sendObject, options);
  }
}

interface Optimiseloss {
  digit: number, 
  tdigit: number ,
  tlen: number, n: number, m: number,
  L: Array<number>,
  U: Array<number>,
  A: Array<number>,
  data: Array<number>,
  gamma: number,
  kappa: number,
  gstrength: number,
  tail: number,
  etLopt: boolean,
  etLmin: number,
  etLmax: number,
  names: Array<string>,
  w: Array<number>,
  initial: Array<number>,
  Q: Array<number>,
  etl: number,
  var: number,
  vaRindex: number,
  breakdown: Array<number>,
  back: number,
  message: string,
  mctr: Array<number>,
  risk: number,
  alpha: Array<number>,
  expreturn: number,
  CVARGLprob: boolean,
  targetReturn: Array<number>,
  loss: number,
  losSopt: boolean,
  losSmin: number,
  losSmax: number,
  result:Result
}

interface Result{
  

  gross: number,
  longvalue: number,
  shortvalue: number,
  shortoverlong: number,
  minhold: number,
  mintrade: number,
  cost: number,
  turnover: number,
  basket: number,
  trades: number,
  var: number,
  vaRindex: number,
  etl: number,
  loss: number,
  breakdown: Array<number>,
  risk: number,
  expreturn: number,
  mctr: Array<number>,
  fmctr: Array<number>,
  sPmctr: Array<number>,
  fx: Array<number>,
  facrisk: number,
  specrisk: number
  
}