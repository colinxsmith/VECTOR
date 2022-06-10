import { Component, Inject, ElementRef } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import * as d3 from 'd3';

@Component({
  selector: 'app-optimiseloss',
  templateUrl: './optimiseloss.component.html',
  styleUrls: ['./optimiseloss.component.css']
})
export class OptimiselossComponent {
  width = 500;
  height = 400;
  format = d3.format('0.6f')
  opt: Array<Optimiseloss> = [];
  constructor(private http: HttpClient, @Inject('BASE_URL') private baseUrl: string, public element: ElementRef) {
    http.get<Optimiseloss[]>(baseUrl + 'optimise/LOSS').subscribe(result => {
      this.opt.push(result[0]);
      console.log(this.opt, result);
    }, error => console.error(error));
  }
  over(e: MouseEvent, inout = false) {
    d3.select(e.target as HTMLInputElement & EventTarget).classed('over', inout);
  }
  sendStep() {
    let back = +(d3.select(this.element.nativeElement).select('input.step.g').node() as HTMLInputElement & Event).value;
    console.log(back, this.opt[0].etLmax);
    this.opt[0].gamma = back;
    back = +(d3.select(this.element.nativeElement).select('input.step.L').node() as HTMLInputElement & Event).value;
    this.opt[0].losSmax = back;
    if (this.opt[0].losSmax != undefined) {
      this.opt[0].losSopt = true;
      this.opt[0].losSmin = back;
    } else {
      this.opt[0].losSopt = false;
    }
    this.sendData('optimise/LOSS', this.opt[0])
      .subscribe(ddd => {
        console.log(ddd);
        this.opt = []
        this.opt.push(ddd[0]);
        console.log(this.opt[0].w);
      }, error => console.error(error));
  }
  sendData(key = 'optimise/LOSS', sendObject: Optimiseloss) {
    const options = {
      headers: new HttpHeaders()
        .set('Content-Type', 'application/json')
    };
    return this.http.post<Array<Optimiseloss>>(`${this.baseUrl}${key}`, sendObject, options);
  }
}

interface Optimiseloss {
  digit: null,
  tdigit: null,
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
  "targetReturn": Array<number>,
  "loss": number,
  "losSopt": boolean,
  "losSmin": number,
  "losSmax": number
}
